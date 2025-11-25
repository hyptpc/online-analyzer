#!/usr/bin/env python3
import time
import re
import sys
import os
import datetime
import sqlite3
import base64
from collections import deque, defaultdict

# +----------------------------+
# | ENVIRONMENT PROXY SETTINGS |
# +----------------------------+
os.environ["NO_PROXY"] = "localhost,127.0.0.1,::1"

# +---------+
# | IMPORTS |
# +---------+
from selenium import webdriver
from selenium.webdriver.firefox.options import Options
from selenium.webdriver.firefox.service import Service
from selenium.webdriver.common.print_page_options import PrintOptions
from webdriver_manager.firefox import GeckoDriverManager
from selenium.webdriver.common.by import By

# +----------------------+
# | LAYOUT CONFIGURATION |
# +----------------------+
MY_LAYOUT = [
    ["RUN",           "Event Number",   "Pi/K"],
    ["----------------", "----------------", "----------------"],
    ["BHT",           "10M-Clock",      "Spill"],
    ["T0",            "TM",             "Real-Time"],
    ["BH2",           "SY",             "Live-Time"],
    ["BAC",           "K-Beam",         "L1-Req"],
    ["KVC",           "Pi-Beam",        "L1-Acc"],
    ["T1",            "Beam",           "L2-Acc"],
    ["HTOF-OR",          "CVC",            "Level1-PS"],
    ["HTOF-NIM-Mp2",      "SAC3",           "TRIG-A-PS"],
    ["HTOF-Mp3",      "SFV",            "TRIG-B-PS"],
    ["HTOF-Fwd",      "Clock-PS",       "TRIG-C-PS"],
    ["BEAM-A",        "TRIG-A",         "TRIG-D-PS"],
    ["BEAM-B",        "TRIG-B",         "TRIG-E-PS"],
    ["BEAM-C",        "TRIG-C",         "TRIG-F-PS"],
    ["BEAM-D",        "TRIG-D",         "TRIG-PSOR-A"],
    ["BEAM-E",        "TRIG-E",         "TRIG-PSOR-B"],
    ["BEAM-F",        "TRIG-F",         "Reserve2-PS"],
    ["----------------", "----------------", "----------------"],
    ["Beam/TM",       "Live/Real",      "DAQ-Eff"],
    ["L1Req/Beam",    "L2-Eff",         "Duty-Factor"],
]

COLUMN_WIDTH = 42
TARGET_URL = "http://kazan.intra.j-parc.jp:8082/?monitoring=50"
TRIGGER_KEY = "10M-Clock"
TRIGGER_THRESHOLD = 20_000_000
CHECK_INTERVAL = 0.1
STABLE_WAIT_SEC = 0.3
DB_FILE = "/home/sks/share/monitor-tmp/scaler_data.db"
  
# Initial Settings
CURRENT_HISTORY_LEN = 10
SETTINGS_FILE = "scaler_settings.txt"

# Initialize History Data
all_history = defaultdict(lambda: deque(maxlen=CURRENT_HISTORY_LEN))
last_trigger_val = -1
last_change_time = time.time()
is_recorded_cycle = False

def clear_screen():
    os.system('cls' if os.name == 'nt' else 'clear')

def check_and_update_settings():
    global CURRENT_HISTORY_LEN, all_history
    if not os.path.exists(SETTINGS_FILE):
        with open(SETTINGS_FILE, "w") as f:
            f.write(str(CURRENT_HISTORY_LEN))
        return
    try:
        with open(SETTINGS_FILE, "r") as f:
            content = f.read().strip()
            if content.isdigit():
                new_len = int(content)
                if new_len > 0 and new_len != CURRENT_HISTORY_LEN:
                    print(f"\n [SETTINGS] Changing average window: {CURRENT_HISTORY_LEN} -> {new_len}")
                    CURRENT_HISTORY_LEN = new_len
                    new_history = defaultdict(lambda: deque(maxlen=CURRENT_HISTORY_LEN))
                    for k, old_deque in all_history.items():
                        new_deque = deque(old_deque, maxlen=CURRENT_HISTORY_LEN)
                        new_history[k] = new_deque
                    all_history = new_history
                    time.sleep(1)
    except Exception as e:
        print(f" [Settings Error] {e}")

def parse_scaler_text(raw_text):
    data = {}
    pattern = r'([A-Za-z0-9\-_/]+(?: [A-Za-z0-9\-_/]+)*)\s*[:=]?\s*([\d,]+\.?\d*|nan|inf|OFF|ON)'
    matches = re.findall(pattern, raw_text)
    for key, val_str in matches:
        key = key.strip()
        try:
            clean_val = val_str.replace(',', '')
            data[key] = float(clean_val)
        except ValueError:
            data[key] = val_str
    return data

def get_display_value(key):
    if key not in all_history or not all_history[key]:
        return "-"
    history = all_history[key]
    sample = history[-1]
    if isinstance(sample, float):
        try:
            valid_nums = [x for x in history if isinstance(x, float)]
            if valid_nums:
                avg_val = sum(valid_nums) / len(history)
                sum_val = sum(valid_nums)
                if avg_val > 100 or avg_val == 0:
                    avg_str = f"{avg_val:,.0f}"
                else:
                    avg_str = f"{avg_val:,.4f}"
                total_str = f"{sum_val:,.1f}"
                return f"{avg_str} ({total_str})"
        except:
            pass
        return str(sample)
    else:
        return str(sample)

# =======================================================
# EXPORT LOGIC
# =======================================================
def prepare_files(trigger_val):
    now_str = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    # 1. Whitelist Creation
    target_keys_set = set()
    for row in MY_LAYOUT:
        for item in row:
            if not item.startswith("---") and not item.startswith("#"):
                target_keys_set.add(item)
    target_keys_set.add(TRIGGER_KEY)
    target_keys = sorted(list(target_keys_set))
    
    row_data = {}
    for k in target_keys:
        if all_history[k]:
            row_data[k] = all_history[k][-1]
        else:
            row_data[k] = None

    # 2. Save to SQLite
    try:
        os.makedirs(os.path.dirname(DB_FILE), exist_ok=True)
        conn = sqlite3.connect(DB_FILE)
        cursor = conn.cursor()
        cursor.execute('CREATE TABLE IF NOT EXISTS scaler_log (timestamp TEXT PRIMARY KEY)')
        
        cursor.execute("PRAGMA table_info(scaler_log)")
        existing_cols = {row[1] for row in cursor.fetchall()}
        
        for k in target_keys:
            if k not in existing_cols:
                try:
                    cursor.execute(f'ALTER TABLE scaler_log ADD COLUMN "{k}" NUMERIC')
                except Exception as ex:
                    print(f"\n[DB Alter Error] {ex}")

        ins_cols = ', '.join(['timestamp'] + [f'"{k}"' for k in target_keys])
        placeholders = ', '.join(['?'] * (1 + len(target_keys)))
        vals = [now_str] + [row_data.get(k) for k in target_keys]
        
        cursor.execute(f'INSERT INTO scaler_log ({ins_cols}) VALUES ({placeholders})', vals)
        conn.commit()
        conn.close()
    except Exception as e:
        print(f"\n[DB Error] {e}")

    # 3. HTML Generation
    print_filename = "scaler_print.html"
    live_filename = "scaler_live.html"
    pdf_filename = "scaler_report_latest.pdf"

    css_style = """
        body { font-family: 'Arial Narrow', sans-serif; padding: 10px; font-size: 11px; }
        .header-container { display: flex; justify-content: space-between; align-items: flex-end; border-bottom: 2px solid #000; margin-bottom: 5px; padding-bottom: 2px; }
        h2 { margin: 0; font-size: 16px; }
        .meta { margin: 0; font-size: 11px; text-align: right; color: #333; }
        table { width: 100%; border-collapse: collapse; border: 1px solid #000; }
        td { border: 1px solid #888; padding: 2px 4px; height: 18px; vertical-align: middle; }
        .label { background-color: #f5f5f5; width: 25%; font-weight: bold; white-space: nowrap; overflow: hidden; }
        .value { text-align: right; font-family: 'Arial Narrow', monospace; font-weight: bold; font-size: 11px; white-space: nowrap; }
        .section-header { background-color: #444; color: #fff; text-align: center; font-weight: bold; font-size: 11px; padding: 2px; }
        .print-btn {
            display: inline-block; padding: 8px 20px; background-color: #007bff; color: white; 
            text-decoration: none; font-weight: bold; border-radius: 4px; margin-bottom: 15px; font-size: 14px;
            box-shadow: 2px 2px 4px rgba(0,0,0,0.2);
        }
        .print-btn:hover { background-color: #0056b3; }
        .no-print { display: none; }
        @page { size: A4 portrait; margin: 5mm; }
    """

    table_html = "<table>"
    for row_items in MY_LAYOUT:
        is_separator = all(item.startswith("---") for item in row_items)
        if is_separator:
             table_html += '<tr style="background-color:#ccc;"><td colspan="6" style="height:2px; padding:0; border:none;"></td></tr>'
             continue
        table_html += "<tr>"
        for key in row_items:
            if key.startswith("---"):
                 table_html += '<td style="border:none; background:#fff;"></td><td style="border:none; background:#fff;"></td>'
            elif key.startswith("#"):
                 table_html += f'<td colspan="2" class="section-header">{key[1:]}</td>'
            else:
                 val_disp = get_display_value(key)
                 table_html += f'<td class="label">{key}</td><td class="value">{val_disp}</td>'
        table_html += "</tr>"
    table_html += "</table>"

    header_html = f"""
        <div class="header-container">
            <h2>ONLINE SCALER</h2>
            <div class="meta">{now_str}<br>Trig: {trigger_val:,.0f} (Avg of {CURRENT_HISTORY_LEN})</div>
        </div>
    """

    print_html = f"<!DOCTYPE html><html><head><meta charset='utf-8'><style>{css_style}</style></head><body>{header_html}{table_html}</body></html>"
    with open(print_filename, "w", encoding="utf-8") as f:
        f.write(print_html)

    live_html = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <meta charset="utf-8">
        <meta http-equiv='refresh' content='3'>
        <title>Live Monitor</title>
        <style>
            {css_style}
            .no-print {{ display: block !important; text-align: center; }}
            @media print {{ .no-print {{ display: none !important; }} }}
        </style>
    </head>
    <body>
        <div class='no-print'>
            <a href='{pdf_filename}' target='_blank' class='print-btn'>🖨️ Open PDF Report</a>
        </div>
        {header_html}
        {table_html}
    </body>
    </html>
    """
    with open(live_filename, "w", encoding="utf-8") as f:
        f.write(live_html)
    
    return os.path.abspath(print_filename)

def generate_pdf(driver, html_path):
    pdf_filename = "scaler_report_latest.pdf"
    original_window = driver.current_window_handle
    try:
        driver.switch_to.new_window('tab')
        driver.get(f"file://{html_path}")
        print_ops = PrintOptions()
        print_ops.orientation = 'portrait'
        print_ops.shrink_to_fit = True
        print_ops.scale = 0.85 
        pdf_b64 = driver.print_page(print_ops)
        with open(pdf_filename, 'wb') as f:
            f.write(base64.b64decode(pdf_b64))
    except Exception as e:
        print(f"\n[PDF Error] {e}")
    finally:
        driver.close()
        driver.switch_to.window(original_window)

def print_user_layout(available_keys):
    clear_screen()
    trig_curr = "Waiting"
    if all_history[TRIGGER_KEY]:
        val = all_history[TRIGGER_KEY][-1]
        trig_curr = f"{val:,.0f}" if isinstance(val, float) else str(val)

    total_width = COLUMN_WIDTH * 3 + 4
    print("=" * total_width)
    print(f"  SCALER DASHBOARD (Avg & Sum) | Trigger: {trig_curr}")
    print("=" * total_width)

    for row_items in MY_LAYOUT:
        line_str = "| "
        for item_key in row_items:
            content = ""
            if item_key.startswith("---"): content = ""
            elif item_key.startswith("#"): content = f"[{item_key[1:]}]"
            else:
                val_disp = get_display_value(item_key)
                content = f"{item_key} : {val_disp}"
            line_str += f"{content:<{COLUMN_WIDTH}} "
        print(line_str)

    print("=" * total_width)
    print(f" [INFO] Reports updated: scaler_live.html, PDF, DB")
    
    layout_keys = set()
    for row in MY_LAYOUT:
        for k in row: layout_keys.add(k)
    missing = [k for k in available_keys if k not in layout_keys]
    if missing:
        print(f"\n[Unused Keys]: {', '.join(missing[:5])} ...")
    print("\n")

def main():
    global last_trigger_val, last_change_time, is_recorded_cycle
    
    # FIREFOX PROXY SETTINGS
    options = Options()
    options.add_argument('--headless')    
    options.set_preference("network.proxy.type", 1) # 1 = Manual
    options.set_preference("network.proxy.http", "buyu.monitor.k18br")
    options.set_preference("network.proxy.http_port", 8080)
    options.set_preference("network.proxy.ssl", "buyu.monitor.k18br")
    options.set_preference("network.proxy.ssl_port", 8080)
    options.set_preference("network.proxy.no_proxies_on", "localhost, 127.0.0.1")

    print(f"Connecting to {TARGET_URL} via Firefox Proxy...")
    driver = webdriver.Firefox(service=Service(GeckoDriverManager().install()), options=options)

    try:
        driver.get(TARGET_URL)
        time.sleep(3)
        check_and_update_settings()
        clear_screen()
        print("Waiting for trigger...")

        while True:
            check_and_update_settings()

            try:
                element = driver.find_element(By.ID, "onlineGUI_drawing")
                raw_text = element.get_attribute("innerText")
                current_data = parse_scaler_text(raw_text)
                
                trigger_raw = current_data.get(TRIGGER_KEY)
                if isinstance(trigger_raw, str): trigger_raw = 0

                if trigger_raw is not None:
                    now = time.time()
                    if trigger_raw == last_trigger_val:
                        elapsed = now - last_change_time
                        
                        if elapsed >= STABLE_WAIT_SEC and not is_recorded_cycle and trigger_raw > TRIGGER_THRESHOLD:
                            
                            # Calc Pi/K
                            pi_val = current_data.get("Pi-Beam", 0)
                            k_val  = current_data.get("K-Beam", 0)
                            if isinstance(pi_val, float) and isinstance(k_val, float) and k_val > 0:
                                current_data["Pi/K"] = pi_val / k_val
                            else:
                                current_data["Pi/K"] = 0.0

                            for k, v in current_data.items():
                                all_history[k].append(v)
                            is_recorded_cycle = True
                            
                            html_path = prepare_files(trigger_raw)
                            generate_pdf(driver, html_path)
                            print_user_layout(current_data.keys())
                    else:
                        last_trigger_val = trigger_raw
                        last_change_time = now
                        is_recorded_cycle = False

                    if not is_recorded_cycle:
                        sys.stdout.write(f"\r >> {TRIGGER_KEY}: {trigger_raw:,.0f} | Stable: {now - last_change_time:.1f}s   ")
                        sys.stdout.flush()
            except Exception:
                pass
            time.sleep(CHECK_INTERVAL)

    except KeyboardInterrupt:
        print("\nStopped.")
    finally:
        driver.quit()

if __name__ == "__main__":
    main()
