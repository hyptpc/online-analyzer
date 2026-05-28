Online-analyzer
============

K1.8BR online analyzer

This explains how to setup and build the online-analyzer from each user’s local PC after connecting to KEKCC.

## Install and build

By default, usr is empty and is created by copying from example.
User part is removed from git management because it will be changed frequently by each user.
First, copy files from example to usr. And make.

```sh
git clone --branch e72 git@github.com:hyptpc/online-analyzer.git
cd online-analyzer
```
Copy the Makefile and setup files:
```sh
./script/copy_makefile.sh
```
Edit src/common.mk and replace:
```sh
unpacker_config := unpacker-config
```
with:
```sh
unpacker_config := /group/had/sks/software/unpacker/e70/bin/unpacker-config
```
Then build the analyzer:
```sh
cd src
make
```

## Parameter setup
The parameter files are managed in a separate repository.
Clone the parameter repository:
```sh
git clone --branch e72 git@github.com:hyptpc/param.git
```
The parameter directory can be placed anywhere depending on each user’s environment.

## Run
Example execution:
```sh
./bin/e72 ../param/conf/analyzer_e72_online.conf /hsm/had/sks/E72/JPARC2025Nov/e72_2025nov/run00001.dat
```
Here:
* `../param/conf/analyzer_e72_online.conf`
    * analyzer configuration file from the param repository
* `/hsm/had/sks/E72/JPARC2025Nov/e72_2025nov/run00001.dat`
    * raw data file path on KEKCC

## Port setting and local host access

The monitor port can be changed by editing the analyzer source code.
For example, in:
```sh
src/analyzer/src/user_raw_e72.cc
```
you may find:
```cpp
int port = 8081;
```

The online monitor will use this port number.
For example:
- `int port = 8081;`
  → access with `http://localhost:8081`
- `int port = 9090;`
  → access with `http://localhost:9090`
  
To access the monitor from your local PC browser, SSH port forwarding is required.
Open another terminal on your local PC and connect to KEKCC with:

```sh
ssh -L 8081:localhost:8081 <your_account>@<server>.cc.kek.jp
```
- Replace `<your_account>` with your KEKCC account name.
- Replace `<server>` with the server you want to use (`cw02`, etc.).

After connecting, open the browser on your local PC:
```sh
http://localhost:8081
```
The port number in:

1. `user_xxx.cc`
2. SSH `-L` option
3. browser `localhost:<port>`

must all match.
