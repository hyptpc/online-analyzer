// -*- C++ -*-

#include "TpcPadHelper.hh"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iterator>
#include <cmath>
#include <stdint.h>
#include <cstdlib>

#include "FuncName.hh"
#include "ConfMan.hh"

namespace
{
const std::string&class_name("TpcPadHelper");
}
const Double_t TpcPadHelper::PadParameter[32][6]=
{ // Layer#, #ofPad, division#, radius#, padLength
  {  0,  48,  14.75,  48, 0,  9.0 },
  {  1,  48,  24.25,  48, 0,  9.0 },
  {  2,  72,  33.75,  72, 0,  9.0 },
  {  3,  96,  43.25,  96, 0,  9.0 },
  {  4, 120,  52.75, 120, 0,  9.0 },
  {  5, 144,  62.25, 144, 0,  9.0 },
  {  6, 168,  71.75, 168, 0,  9.0 },
  {  7, 192,  81.25, 192, 0,  9.0 },
  {  8, 216,  90.75, 216, 0,  9.0 },
  {  9, 240, 100.25, 240, 0,  9.0 },
  { 10, 208, 111.50, 241, 0, 12.5 },
  { 11, 218, 124.50, 271, 0, 12.5 },
  { 12, 230, 137.50, 300, 0, 12.5 },
  { 13, 214, 150.50, 330, 0, 12.5 },
  { 14, 212, 163.50, 360, 0, 12.5 },
  { 15, 214, 176.50, 390, 0, 12.5 },
  { 16, 220, 189.50, 420, 0, 12.5 },
  { 17, 224, 202.50, 449, 0, 12.5 },
  { 18, 232, 215.50, 479, 0, 12.5 },
  { 19, 238, 228.50, 509, 0, 12.5 },
  { 20, 244, 241.50, 539, 0, 12.5 },
  { 21, 232, 254.50, 569, 0, 12.5 },
  { 22, 218, 267.50, 599, 0, 12.5 },
  { 23, 210, 280.50, 628, 0, 12.5 },
  { 24, 206, 293.50, 658, 0, 12.5 },
  { 25, 202, 306.50, 688, 0, 12.5 },
  { 26, 200, 319.50, 718, 0, 12.5 },
  { 27, 196, 332.50, 748, 0, 12.5 },
  { 28, 178, 345.50, 777, 0, 12.5 },
  { 29, 130, 358.50, 807, 0, 12.5 },
  { 30, 108, 371.50, 837, 0, 12.5 },
  { 31,  90, 384.50, 867, 0, 12.5 } };

//_____________________________________________________________________________
const Int_t TpcPadHelper::padOnCenterFrame[633] =
{

//Pads on the frame
964,965,966,967,968,969,970,971,972,973,1018,1019,1020,1021,1022,1023,1024,1025,1026,1027,1176,1177,1178,1179,1180,1181,1182,1183,1184,1185,1186,1187,1188,1189,1190,1191,1192,1193,1194,1195,1196,1197,1198,1199,1200,1201,1202,1203,1204,1205,1206,1207,1208,1209,1210,1211,1212,1236,1237,1238,1239,1240,1241,1242,1243,1244,1245,1246,1247,1248,1249,1250,1251,1252,1253,1254,1255,1256,1257,1258,1259,1260,1261,1262,1263,1264,1265,1266,1267,1268,1269,1270,1271,1394,1395,1396,1397,1398,1399,1400,1404,1405,1406,1407,1408,1409,1410,1411,1412,1413,1414,1415,1416,1417,1418,1419,1420,1421,1422,1423,1424,1425,1426,1427,1428,1429,1430,1431,1435,1436,1437,1438,1439,1440,1441,1454,1455,1456,1457,1458,1459,1460,1464,1465,1466,1467,1468,1469,1470,1471,1472,1473,1474,1475,1476,1477,1478,1479,1480,1481,1482,1483,1484,1485,1486,1487,1488,1489,1490,1491,1495,1496,1497,1498,1499,1500,1501,1594,1595,1596,1597,1603,1604,1605,1606,1607,1646,1647,1648,1649,1650,1651,1657,1658,1659,1662,1663,1664,1670,1671,1672,1673,1674,1675,1714,1715,1716,1717,1718,1724,1725,1726,1727,1806,1807,1808,1809,1810,1811,1812,1813,1814,1815,1816,1817,1878,1879,1880,1881,1888,1889,1890,1891,1952,1953,1954,1955,1956,1957,1958,1959,1960,1961,1962,1963,2017,2018,2019,2024,2025,2026,2027,2099,2100,2101,2102,2103,2110,2111,2112,2113,2114,2186,2187,2188,2189,2194,2195,2196,2218,2219,2220,2221,2222,2223,2224,2225,2226,2227,2308,2309,2310,2316,2317,2322,2323,2329,2330,2331,2412,2413,2414,2415,2416,2417,2418,2419,2420,2421,2426,2427,2428,2429,2517,2518,2519,2520,2521,2522,2523,2524,2525,2526,2539,2540,2541,2542,2543,2544,2545,2546,2547,2548,2636,2637,2638,2639,2731,2732,2738,2739,2760,2761,2767,2768,2949,2950,2951,2952,2953,2954,2955,2956,2957,2986,2987,2988,2989,2990,2991,2992,2993,2994,3173,3174,3175,3176,3177,3178,3179,3180,3181,3218,3219,3220,3221,3222,3223,3224,3225,3226,3404,3405,3406,3407,3408,3409,3410,3411,3412,3457,3458,3459,3460,3461,3462,3463,3464,3465,3641,3642,3643,3644,3645,3646,3647,3648,3649,3702,3703,3704,3705,3706,3707,3708,3709,3710,3876,3877,3878,3879,3880,3881,3882,3883,3944,3945,3946,3947,3948,3949,3950,3951,4097,4098,4104,4173,4179,4180,4307,4308,4309,4310,4311,4312,4313,4314,4315,4390,4391,4392,4393,4394,4395,4396,4397,4398,4511,4512,4519,4602,4609,4610,4712,4713,4714,4715,4716,4717,4718,4719,4810,4811,4812,4813,4814,4815,4816,4817,4909,4910,4911,4912,4913,4914,4915,4916,5015,5016,5017,5018,5019,5020,5021,5022,5103,5104,5107,5110,5111,5216,5217,5220,5223,5224,5286,5287,5288,5289,5290,5291,5292,5293,5294,5407,5408,5409,5410,5411,5412,5413,5414,5415,5440,5441,5442,5443,5444,5565,5566,5567,5568,5569,

//Empty Pads
1016,1017,1393,1401,1402,1403,1432,1433,1434,1461,1462,1463,1492,1493,1494,1502,1598,1599,1600,1601,1602,1652,1653,1654,1655,1656,1665,1666,1667,1668,1669,1719,1720,1721,1722,1723,1876,1877,1882,1883,1884,1885,1886,1887,1892,1893,2020,2021,2022,2023,2104,2105,2106,2107,2108,2109,2190,2191,2192,2193,2311,2312,2313,2314,2315,2324,2325,2326,2327,2328,2733,2734,2735,2736,2737,2762,2763,2764,2765,2766,4099,4100,4101,4102,4103,4174,4175,4176,4177,4178,4513,4514,4515,4516,4517,4518,4603,4604,4605,4606,4607,4608,5105,5106,5108,5109,5218,5219,5221,5222

};

//_____________________________________________________________________________
TpcPadParam::TpcPadParam( Int_t asad, Int_t aget, Int_t channel,
                          Int_t pad_id, Int_t layer_id, Int_t row_id )
  : m_asad( asad ),
    m_aget( aget ),
    m_channel( channel ),
    m_pad_id( pad_id ),
    m_layer_id( layer_id ),
    m_row_id( row_id ),
    m_key( Key( asad, aget, channel ) )
{
}

//_____________________________________________________________________________
TpcPadParam::~TpcPadParam( void )
{
}

//_____________________________________________________________________________
void
TpcPadParam::Print( void ) const
{
  std::cout << " AsAd=" << m_asad
	    << " AGet=" << m_aget
	    << " Channel=" << m_channel
	    << " PadId=" << m_pad_id
	    << " LayerId=" << m_layer_id
	    << " RowId=" << m_row_id
	    << " Key=" << m_key
	    << std::endl;
}


//_____________________________________________________________________________
TpcPadHelper::TpcPadHelper( void )
  : m_map(), m_file_name("")
{
}

//_____________________________________________________________________________
TpcPadHelper::~TpcPadHelper( void )
{
}

//_____________________________________________________________________________
Int_t
TpcPadHelper::GetPadId( Int_t asad, Int_t aget, Int_t channel ) const
{
  auto param = GetParam( asad, aget, channel );
  if( param && param->AsAdId()>=0 )
    return param->PadId();
  else
    return -1;
}

//0-origin
//_____________________________________________________________________________
Int_t
TpcPadHelper::GetPadId(Int_t layer_id, Int_t row_id) const
{
 
  if (row_id < 0) return -1;

  Int_t pad_id = 0;
  for (Int_t layer = 0; layer < layer_id; layer++)
    pad_id += static_cast<Int_t>(PadParameter[layer][kNumOfPad]);
  pad_id += row_id;
  return pad_id;
}

//_____________________________________________________________________________
Int_t
TpcPadHelper::GetASADId(Int_t layer, Int_t row) const
{
  Int_t flag = layer/4;
  Int_t section;
  if(flag==0) section=0; //layer 0~3
  else if(flag==1) section=1; //layer 4~7
  else if(layer==30||layer==31) section=3; //layer 30~31
  else section=2; //layer 8~29

  Int_t half = PadParameter[layer][1]/2;
  Int_t division1 = PadParameter[layer][1]/6;
  Int_t division2 = PadParameter[layer][1]*5/6;

  switch(section){
  case 0:
    if(row<half)
      return 0;
    else
      return 1;
  case 1:
    Int_t dummy;
    if(layer%4<2) dummy=2;
    if(2<=layer%4) dummy=5;
    if(row<division1||division2<=row)
      return dummy;
    else if(division1<=row&&row<half)
      return dummy+1;
    else
      return dummy+2;
  case 3:
    return 30;
  default:
    if(row<half)
      return layer-layer%2;
    else
      return layer-layer%2+1;
  }
}

//_____________________________________________________________________________
Int_t
TpcPadHelper::GetAGETId(Int_t asad, Int_t layer, Int_t row) const
{
#ifdef PAD_HELPER_DEBUG
  ValidateRow(layer, row, __func__);
#endif
  Int_t flag=-1;
  switch(asad){
  case 0:
    //AGET-0
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==1 && row==22-ch*2)      flag=0;
    for(Int_t ch=13; ch<=21; ++ch) if(layer==3 && row==46-(ch-13)*2) flag=0;
    for(Int_t ch=23; ch<=37; ++ch) if(layer==3 && row==28-(ch-23)*2) flag=0;
    if(layer==1 && row==0) flag=0;
    //AGET-1
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==1 && row==23-ch*2)      flag=1;
    for(Int_t ch=13; ch<=21; ++ch) if(layer==3 && row==47-(ch-13)*2) flag=1;
    for(Int_t ch=23; ch<=37; ++ch) if(layer==3 && row==29-(ch-23)*2) flag=1;
    if(layer==1 && row==1) flag=1;
    //AGET-2
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==0 && row==22-ch*2)      flag=2;
    for(Int_t ch=13; ch<=21; ++ch) if(layer==2 && row==34-(ch-13)*2) flag=2;
    for(Int_t ch=23; ch<=31; ++ch) if(layer==2 && row==16-(ch-23)*2) flag=2;
    if(layer==0 && row==0) flag=2;
    //AGET-3
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==0 && row==23-ch*2)      flag=3;
    for(Int_t ch=13; ch<=21; ++ch) if(layer==2 && row==35-(ch-13)*2) flag=3;
    for(Int_t ch=23; ch<=31; ++ch) if(layer==2 && row==17-(ch-23)*2) flag=3;
    if(layer==0 && row==1) flag=3;
    return flag;
  case 1:
    //AGET-0
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==1 && row==47-ch*2)      flag=0;
    for(Int_t ch=13; ch<=21; ++ch) if(layer==3 && row==95-(ch-13)*2) flag=0;
    for(Int_t ch=23; ch<=37; ++ch) if(layer==3 && row==77-(ch-23)*2) flag=0;
    if(layer==1 && row==25) flag=0;
    //AGET-1
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==1 && row==46-ch*2)      flag=1;
    for(Int_t ch=13; ch<=21; ++ch) if(layer==3 && row==94-(ch-13)*2) flag=1;
    for(Int_t ch=23; ch<=37; ++ch) if(layer==3 && row==76-(ch-23)*2) flag=1;
    if(layer==1 && row==24) flag=1;
    //AGET-2
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==0 && row==47-ch*2)      flag=2;
    for(Int_t ch=13; ch<=21; ++ch) if(layer==2 && row==71-(ch-13)*2) flag=2;
    for(Int_t ch=23; ch<=37; ++ch) if(layer==2 && row==53-(ch-23)*2) flag=2;
    if(layer==0 && row==25) flag=2;
    //AGET-3
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==0 && row==46-ch*2)      flag=3;
    for(Int_t ch=13; ch<=21; ++ch) if(layer==2 && row==70-(ch-13)*2) flag=3;
    for(Int_t ch=23; ch<=31; ++ch) if(layer==2 && row==52-(ch-23)*2) flag=3;
    if(layer==0 && row==24) flag=3;
    return flag;
  case 2:
    //AGET-0
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==5 && row==22-ch*2)       flag=0;
    for(Int_t ch=13; ch<=21; ++ch) if(layer==5 && row==143-(ch-13)*2) flag=0;
    for(Int_t ch=23; ch<=25; ++ch) if(layer==5 && row==125-(ch-23)*2) flag=0;
    if(layer==5 && row==0) flag=0;
    //AGET-1
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==5 && row==23-ch*2)       flag=1;
    for(Int_t ch=13; ch<=21; ++ch) if(layer==5 && row==142-(ch-13)*2) flag=1;
    for(Int_t ch=23; ch<=25; ++ch) if(layer==5 && row==124-(ch-23)*2) flag=1;
    if(layer==5 && row==1) flag=1;
    //AGET-2
    for(Int_t ch= 0; ch<= 9; ++ch) if(layer==4 && row==18-ch*2)       flag=2;
    for(Int_t ch=12; ch<=20; ++ch) if(layer==4 && row==117-(ch-12)*2) flag=2;
    if(layer==4 && row==119) flag=2;
    //AGET-3
    for(Int_t ch= 0; ch<= 9; ++ch) if(layer==4 && row==19-ch*2)       flag=3;
    for(Int_t ch=12; ch<=20; ++ch) if(layer==4 && row==116-(ch-12)*2) flag=3;
    if(layer==4 && row==118) flag=3;
    return flag;
  case 3:
    //AGET-0
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==5 && row==24+ch*2)      flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==5 && row==46+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=25; ++ch) if(layer==5 && row==66+(ch-23)*2) flag=0;
    //AGET-1
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==5 && row==25+ch*2)      flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==5 && row==47+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=25; ++ch) if(layer==5 && row==67+(ch-23)*2) flag=1;
    //AGET-2
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==4 && row==20+ch*2)      flag=2;
    for(Int_t ch=12; ch<=20; ++ch) if(layer==4 && row==42+(ch-12)*2) flag=2;
    //AGET-3
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==4 && row==21+ch*2)      flag=3;
    for(Int_t ch=12; ch<=20; ++ch) if(layer==4 && row==43+(ch-12)*2) flag=3;
    return flag;
  case 4:
    //AGET-0
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==5 && row==73+ch*2)       flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==5 && row==95+(ch-12)*2)  flag=0;
    for(Int_t ch=23; ch<=25; ++ch) if(layer==5 && row==115+(ch-23)*2) flag=0;
    //AGET-1
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==5 && row==72+ch*2)       flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==5 && row==94+(ch-12)*2)  flag=1;
    for(Int_t ch=23; ch<=25; ++ch) if(layer==5 && row==114+(ch-23)*2) flag=1;
    //AGET-2
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==4 && row==61+ch*2)      flag=2;
    for(Int_t ch=12; ch<=20; ++ch) if(layer==4 && row==83+(ch-12)*2) flag=2;
    //AGET-3
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==4 && row==60+ch*2)      flag=3;
    for(Int_t ch=12; ch<=20; ++ch) if(layer==4 && row==82+(ch-12)*2) flag=3;
    return flag;
  case 5:
    //AGET-0
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==7 && row==30-ch*2)       flag=0;
    for(Int_t ch=12; ch<=16; ++ch) if(layer==7 && row==8-(ch-12)*2)   flag=0;
    for(Int_t ch=17; ch<=21; ++ch) if(layer==7 && row==191-(ch-17)*2) flag=0;
    for(Int_t ch=23; ch<=33; ++ch) if(layer==7 && row==181-(ch-23)*2) flag=0;
    //AGET-1
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==7 && row==31-ch*2)       flag=1;
    for(Int_t ch=12; ch<=16; ++ch) if(layer==7 && row==9-(ch-12)*2)   flag=1;
    for(Int_t ch=17; ch<=21; ++ch) if(layer==7 && row==190-(ch-17)*2) flag=1;
    for(Int_t ch=23; ch<=33; ++ch) if(layer==7 && row==180-(ch-23)*2) flag=1;
    //AGET-2
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==6 && row==26-ch*2)       flag=2;
    for(Int_t ch=12; ch<=14; ++ch) if(layer==6 && row==4-(ch-12)*2)   flag=2;
    for(Int_t ch=15; ch<=21; ++ch) if(layer==6 && row==167-(ch-15)*2) flag=2;
    for(Int_t ch=23; ch<=29; ++ch) if(layer==6 && row==153-(ch-23)*2) flag=2;
    //AGET-3
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==6 && row==27-ch*2)       flag=3;
    for(Int_t ch=12; ch<=14; ++ch) if(layer==6 && row==5-(ch-12)*2)   flag=3;
    for(Int_t ch=15; ch<=21; ++ch) if(layer==6 && row==166-(ch-15)*2) flag=3;
    for(Int_t ch=23; ch<=29; ++ch) if(layer==6 && row==152-(ch-23)*2) flag=3;
    return flag;
  case 6:
    //AGET-0
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==7 && row==32+ch*2)      flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==7 && row==54+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=33; ++ch) if(layer==7 && row==74+(ch-23)*2) flag=0;
    //AGET-1
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==7 && row==33+ch*2)      flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==7 && row==55+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=33; ++ch) if(layer==7 && row==75+(ch-23)*2) flag=1;
    //AGET-2
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==6 && row==28+ch*2)      flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==6 && row==50+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=29; ++ch) if(layer==6 && row==70+(ch-23)*2) flag=2;
    //AGET-3
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==6 && row==29+ch*2)      flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==6 && row==51+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=29; ++ch) if(layer==6 && row==71+(ch-23)*2) flag=3;
    return flag;
  case 7:
    //AGET-0
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==7 && row==97+ch*2)       flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==7 && row==119+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=33; ++ch) if(layer==7 && row==139+(ch-23)*2) flag=0;
    //AGET-1
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==7 && row==96+ch*2)       flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==7 && row==118+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=33; ++ch) if(layer==7 && row==138+(ch-23)*2) flag=1;
    //AGET-2
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==6 && row==85+ch*2)       flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==6 && row==107+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=29; ++ch) if(layer==6 && row==127+(ch-23)*2) flag=2;
    //AGET-3
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==6 && row==84+ch*2)       flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==6 && row==106+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=29; ++ch) if(layer==6 && row==126+(ch-23)*2) flag=3;
    return flag;
  case 8:
    //AGET-0
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==9 && row==ch*2)          flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==9 && row==22+(ch-12)*2)  flag=0;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==9 && row==42+(ch-23)*2)  flag=0;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==9 && row==86+(ch-46)*2)  flag=0;
    for(Int_t ch=57; ch<=63; ++ch) if(layer==9 && row==106+(ch-57)*2) flag=0;
    //AGET-1
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==9 && row==1+ch*2)        flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==9 && row==23+(ch-12)*2)  flag=1;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==9 && row==43+(ch-23)*2)  flag=1;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==9 && row==87+(ch-46)*2)  flag=1;
    for(Int_t ch=57; ch<=63; ++ch) if(layer==9 && row==107+(ch-57)*2) flag=1;
    //AGET-2
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==8 && row==ch*2)         flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==8 && row==22+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==8 && row==42+(ch-23)*2) flag=2;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==8 && row==86+(ch-46)*2) flag=2;
    if(layer==8 && row==106) flag=2;
    //AGET-3
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==8 && row==1+ch*2)       flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==8 && row==23+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==8 && row==43+(ch-23)*2) flag=3;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==8 && row==87+(ch-46)*2) flag=3;
    if(layer==8 && row==107) flag=3;
    return flag;
  case 9:
    //AGET-0
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==9 && row==121+ch*2)      flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==9 && row==143+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==9 && row==163+(ch-23)*2) flag=0;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==9 && row==207+(ch-46)*2) flag=0;
    for(Int_t ch=57; ch<=63; ++ch) if(layer==9 && row==227+(ch-57)*2) flag=0;
    //AGET-1
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==9 && row==120+ch*2)      flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==9 && row==142+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==9 && row==162+(ch-23)*2) flag=1;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==9 && row==206+(ch-46)*2) flag=1;
    for(Int_t ch=57; ch<=63; ++ch) if(layer==9 && row==226+(ch-57)*2) flag=1;
    //AGET-2
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==8 && row==109+ch*2)      flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==8 && row==131+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==8 && row==151+(ch-23)*2) flag=2;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==8 && row==195+(ch-46)*2) flag=2;
    if(layer==8 && row==215) flag=2;
    //AGET-3
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==8 && row==108+ch*2)      flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==8 && row==130+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==8 && row==150+(ch-23)*2) flag=3;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==8 && row==194+(ch-46)*2) flag=3;
    if(layer==8 && row==214) flag=3;
    return flag;
  case 10:
    //AGET-0
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==11 && row==1+ch*2)       flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==11 && row==23+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==11 && row==43+(ch-23)*2) flag=0;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==11 && row==87+(ch-46)*2) flag=0;
    if(layer==11 && row==107) flag=0;
    //AGET-1
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==11 && row==ch*2)          flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==11 && row==22+(ch-12)*2)  flag=1;
    for(Int_t ch=23; ch<=24; ++ch) if(layer==11 && row==42+(ch-23)*2)  flag=1;
    for(Int_t ch=28; ch<=44; ++ch) if(layer==11 && row==52+(ch-28)*2)  flag=1;
    for(Int_t ch=46; ch<=52; ++ch) if(layer==11 && row==86+(ch-46)*2)  flag=1;
    for(Int_t ch=57; ch<=58; ++ch) if(layer==11 && row==106+(ch-57)*2) flag=1;
    //AGET-2
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==10 && row==ch*2)         flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==10 && row==22+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=30; ++ch) if(layer==10 && row==42+(ch-23)*2) flag=2;
    for(Int_t ch=32; ch<=44; ++ch) if(layer==10 && row==60+(ch-32)*2) flag=2;
    for(Int_t ch=49; ch<=54; ++ch) if(layer==10 && row==92+(ch-49)*2) flag=2;
    if(layer==10 && row==86) flag=2;
    //AGET-3
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==10 && row==1+ch*2)       flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==10 && row==23+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=25; ++ch) if(layer==10 && row==43+(ch-23)*2) flag=3;
    for(Int_t ch=27; ch<=29; ++ch) if(layer==10 && row==51+(ch-27)*2) flag=3;
    for(Int_t ch=32; ch<=44; ++ch) if(layer==10 && row==61+(ch-32)*2) flag=3;
    for(Int_t ch=48; ch<=54; ++ch) if(layer==10 && row==91+(ch-48)*2) flag=3;
    if(layer==10 && row==87) flag=3;
    return flag;
  case 11:
    //AGET-0
    for(Int_t ch= 0; ch<= 1; ++ch) if(layer==11 && row==110+ch*2)      flag=0;
    for(Int_t ch= 4; ch<=10; ++ch) if(layer==11 && row==118+(ch-4)*2)  flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==11 && row==132+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=30; ++ch) if(layer==11 && row==152+(ch-23)*2) flag=0;
    for(Int_t ch=33; ch<=44; ++ch) if(layer==11 && row==172+(ch-33)*2) flag=0;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==11 && row==196+(ch-46)*2) flag=0;
    if(layer==11 && row==216) flag=0;
    //AGET-1
    for(Int_t ch= 0; ch<= 1; ++ch) if(layer==11 && row==109+ch*2)      flag=1;
    for(Int_t ch= 5; ch<=10; ++ch) if(layer==11 && row==119+(ch-5)*2)  flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==11 && row==131+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=30; ++ch) if(layer==11 && row==151+(ch-23)*2) flag=1;
    for(Int_t ch=34; ch<=44; ++ch) if(layer==11 && row==173+(ch-34)*2) flag=1;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==11 && row==195+(ch-46)*2) flag=1;
    for(Int_t ch=57; ch<=58; ++ch) if(layer==11 && row==215+(ch-57)*2) flag=1;
    //AGET-2
    for(Int_t ch= 0; ch<= 5; ++ch) if(layer==10 && row==105+ch*2)      flag=2;
    for(Int_t ch= 8; ch<=10; ++ch) if(layer==10 && row==121+(ch-8)*2)  flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==10 && row==127+(ch-12)*2) flag=2;
    for(Int_t ch=25; ch<=44; ++ch) if(layer==10 && row==151+(ch-25)*2) flag=2;
    for(Int_t ch=46; ch<=54; ++ch) if(layer==10 && row==191+(ch-46)*2) flag=2;
    if(layer==10 && row==147) flag=2;
    //AGET-3
    for(Int_t ch= 0; ch<= 6; ++ch) if(layer==10 && row==104+ch*2)      flag=3;
    for(Int_t ch= 8; ch<=10; ++ch) if(layer==10 && row==120+(ch-8)*2)  flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==10 && row==126+(ch-12)*2) flag=3;
    for(Int_t ch=26; ch<=28; ++ch) if(layer==10 && row==152+(ch-26)*2) flag=3;
    for(Int_t ch=30; ch<=44; ++ch) if(layer==10 && row==160+(ch-30)*2) flag=3;
    for(Int_t ch=46; ch<=54; ++ch) if(layer==10 && row==190+(ch-46)*2) flag=3;
    if(layer==10 && row==146) flag=3;
    return flag;
  case 12:
    //AGET-0
    for(Int_t ch= 0; ch<= 9; ++ch) if(layer==13 && row==1+ch*2)       flag=0;
    for(Int_t ch=13; ch<=21; ++ch) if(layer==13 && row==25+(ch-13)*2) flag=0;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==13 && row==43+(ch-23)*2) flag=0;
    for(Int_t ch=46; ch<=54; ++ch) if(layer==13 && row==87+(ch-46)*2) flag=0;
    //AGET-1
    for(Int_t ch= 0; ch<= 9; ++ch) if(layer==13 && row==ch*2)         flag=1;
    for(Int_t ch=13; ch<=21; ++ch) if(layer==13 && row==24+(ch-13)*2) flag=1;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==13 && row==42+(ch-23)*2) flag=1;
    for(Int_t ch=46; ch<=54; ++ch) if(layer==13 && row==86+(ch-46)*2) flag=1;
    //AGET-2
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==12 && row==1+ch*2)        flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==12 && row==23+(ch-12)*2)  flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==12 && row==43+(ch-23)*2)  flag=2;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==12 && row==87+(ch-46)*2)  flag=2;
    for(Int_t ch=58; ch<=59; ++ch) if(layer==12 && row==109+(ch-58)*2) flag=2;
    //AGET-3
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==12 && row==ch*2)          flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==12 && row==22+(ch-12)*2)  flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==12 && row==42+(ch-23)*2)  flag=3;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==12 && row==86+(ch-46)*2)  flag=3;
    for(Int_t ch=58; ch<=59; ++ch) if(layer==12 && row==108+(ch-58)*2) flag=3;
    return flag;
  case 13:
    //AGET-0
    for(Int_t ch= 1; ch<=10; ++ch) if(layer==13 && row==110+(ch-1)*2)  flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==13 && row==130+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=42; ++ch) if(layer==13 && row==150+(ch-23)*2) flag=0;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==13 && row==194+(ch-46)*2) flag=0;
    //AGET-1
    for(Int_t ch= 2; ch<=10; ++ch) if(layer==13 && row==111+(ch-2)*2)  flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==13 && row==129+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=43; ++ch) if(layer==13 && row==149+(ch-23)*2) flag=1;
    for(Int_t ch=47; ch<=55; ++ch) if(layer==13 && row==195+(ch-47)*2) flag=1;
    if(layer==13 && row==213) flag=1;
    //AGET-2
    for(Int_t ch= 4; ch<=10; ++ch) if(layer==12 && row==124+(ch-4)*2)  flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==12 && row==138+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==12 && row==158+(ch-23)*2) flag=2;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==12 && row==202+(ch-46)*2) flag=2;
    for(Int_t ch=57; ch<=60; ++ch) if(layer==12 && row==222+(ch-57)*2) flag=2;
    if(layer==12 && row==118) flag=2;
    if(layer==12 && row==120) flag=2;
    //AGET-3
    for(Int_t ch= 5; ch<=10; ++ch) if(layer==12 && row==125+(ch-5)*2)  flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==12 && row==137+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==12 && row==157+(ch-23)*2) flag=3;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==12 && row==201+(ch-46)*2) flag=3;
    for(Int_t ch=57; ch<=61; ++ch) if(layer==12 && row==221+(ch-57)*2) flag=3;
    if(layer==12 && row==119) flag=3;
    if(layer==12 && row==121) flag=3;
    return flag;
  case 14:
    //AGET-0
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==15 && row==1+ch*2)       flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==15 && row==23+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==15 && row==43+(ch-23)*2) flag=0;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==15 && row==87+(ch-46)*2) flag=0;
    //AGET-1
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==15 && row==ch*2)         flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==15 && row==22+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==15 && row==42+(ch-23)*2) flag=1;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==15 && row==86+(ch-46)*2) flag=1;
    if(layer==15 && row==106) flag=1;
    //AGET-2
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==14 && row==ch*2)          flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==14 && row==22+(ch-12)*2)  flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==14 && row==42+(ch-23)*2)  flag=2;
    for(Int_t ch=46; ch<=51; ++ch) if(layer==14 && row==86+(ch-46)*2)  flag=2;
    for(Int_t ch=54; ch<=55; ++ch) if(layer==14 && row==102+(ch-54)*2) flag=2;
    //AGET-3
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==14 && row==1+ch*2)        flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==14 && row==23+(ch-12)*2)  flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==14 && row==43+(ch-23)*2)  flag=3;
    for(Int_t ch=46; ch<=50; ++ch) if(layer==14 && row==87+(ch-46)*2)  flag=3;
    for(Int_t ch=54; ch<=55; ++ch) if(layer==14 && row==103+(ch-54)*2) flag=3;
    return flag;
  case 15:
    //AGET-0
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==15 && row==108+ch*2)      flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==15 && row==130+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==15 && row==150+(ch-23)*2) flag=0;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==15 && row==194+(ch-46)*2) flag=0;
    //AGET-1
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==15 && row==107+ch*2)      flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==15 && row==129+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==15 && row==149+(ch-23)*2) flag=1;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==15 && row==193+(ch-46)*2) flag=1;
    if(layer==15 && row==213) flag=1;
    //AGET-2
    if(layer==14 && row==107) flag=2;
    if(layer==14 && row==109) flag=2;
    for(Int_t ch= 4; ch<=10; ++ch) if(layer==14 && row==115+(ch-4)*2)  flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==14 && row==129+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==14 && row==149+(ch-23)*2) flag=2;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==14 && row==193+(ch-46)*2) flag=2;
    //AGET-3
    if(layer==14 && row==106) flag=3;
    if(layer==14 && row==108) flag=3;
    for(Int_t ch= 5; ch<=10; ++ch) if(layer==14 && row==116+(ch-5)*2)  flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==14 && row==128+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==14 && row==148+(ch-23)*2) flag=3;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==14 && row==192+(ch-46)*2) flag=3;
    return flag;
  case 16:
    //AGET-0
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==17 && row==ch*2)          flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==17 && row==22+(ch-12)*2)  flag=0;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==17 && row==42+(ch-23)*2)  flag=0;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==17 && row==86+(ch-46)*2)  flag=0;
    for(Int_t ch=57; ch<=59; ++ch) if(layer==17 && row==106+(ch-57)*2) flag=0;
    //AGET-1
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==17 && row==1+ch*2)        flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==17 && row==23+(ch-12)*2)  flag=1;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==17 && row==43+(ch-23)*2)  flag=1;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==17 && row==87+(ch-46)*2)  flag=1;
    for(Int_t ch=57; ch<=59; ++ch) if(layer==17 && row==107+(ch-57)*2) flag=1;
    //AGET-2
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==16 && row==ch*2)          flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==16 && row==22+(ch-12)*2)  flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==16 && row==42+(ch-23)*2)  flag=2;
    for(Int_t ch=46; ch<=49; ++ch) if(layer==16 && row==86+(ch-46)*2)  flag=2;
    for(Int_t ch=52; ch<=55; ++ch) if(layer==16 && row==98+(ch-52)*2)  flag=2;
    for(Int_t ch=57; ch<=58; ++ch) if(layer==16 && row==106+(ch-57)*2) flag=2;
    //AGET-3
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==16 && row==1+ch*2)       flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==16 && row==23+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==16 && row==43+(ch-23)*2) flag=3;
    for(Int_t ch=46; ch<=48; ++ch) if(layer==16 && row==87+(ch-46)*2) flag=3;
    for(Int_t ch=52; ch<=55; ++ch) if(layer==16 && row==99+(ch-52)*2) flag=3;
    if(layer==16 && row== 10) flag=3;
    if(layer==16 && row==109) flag=3;
    return flag;
  case 17:
    //AGET-0
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==17 && row==113+ch*2)      flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==17 && row==135+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==17 && row==155+(ch-23)*2) flag=0;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==17 && row==199+(ch-46)*2) flag=0;
    for(Int_t ch=57; ch<=59; ++ch) if(layer==17 && row==219+(ch-57)*2) flag=0;
    //AGET-1
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==17 && row==112+ch*2)      flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==17 && row==134+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==17 && row==154+(ch-23)*2) flag=1;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==17 && row==198+(ch-46)*2) flag=1;
    for(Int_t ch=57; ch<=59; ++ch) if(layer==17 && row==218+(ch-57)*2) flag=1;
    //AGET-2
    for(Int_t ch= 0; ch<= 5; ++ch) if(layer==16 && row==111+ch*2)      flag=2;
    for(Int_t ch= 8; ch<=10; ++ch) if(layer==16 && row==127+(ch-8)*2)  flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==16 && row==133+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==16 && row==153+(ch-23)*2) flag=2;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==16 && row==197+(ch-46)*2) flag=2;
    if(layer==16 && row==217) flag=2;
    if(layer==16 && row==219) flag=2;
    //AGET-3
    for(Int_t ch= 0; ch<= 5; ++ch) if(layer==16 && row==110+ch*2)      flag=3;
    for(Int_t ch= 9; ch<=10; ++ch) if(layer==16 && row==128+(ch-9)*2)  flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==16 && row==132+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==16 && row==152+(ch-23)*2) flag=3;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==16 && row==196+(ch-46)*2) flag=3;
    for(Int_t ch=57; ch<=58; ++ch) if(layer==16 && row==216+(ch-57)*2) flag=3;
    return flag;
  case 18:
    //AGET-0
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==19 && row==1+ch*2)        flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==19 && row==23+(ch-12)*2)  flag=0;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==19 && row==43+(ch-23)*2)  flag=0;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==19 && row==87+(ch-46)*2)  flag=0;
    for(Int_t ch=57; ch<=62; ++ch) if(layer==19 && row==107+(ch-57)*2) flag=0;
    //AGET-1
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==19 && row==ch*2)          flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==19 && row==22+(ch-12)*2)  flag=1;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==19 && row==42+(ch-23)*2)  flag=1;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==19 && row==86+(ch-46)*2)  flag=1;
    for(Int_t ch=57; ch<=63; ++ch) if(layer==19 && row==106+(ch-57)*2) flag=1;
    //AGET-2
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==18 && row==ch*2)          flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==18 && row==22+(ch-12)*2)  flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==18 && row==42+(ch-23)*2)  flag=2;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==18 && row==86+(ch-46)*2)  flag=2;
    for(Int_t ch=57; ch<=61; ++ch) if(layer==18 && row==106+(ch-57)*2) flag=2;
    //AGET-3
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==18 && row==1+ch*2)        flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==18 && row==23+(ch-12)*2)  flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==18 && row==43+(ch-23)*2)  flag=3;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==18 && row==87+(ch-46)*2)  flag=3;
    for(Int_t ch=57; ch<=61; ++ch) if(layer==18 && row==107+(ch-57)*2) flag=3;
    return flag;
  case 19:
    //AGET-0
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==19 && row==120+ch*2)      flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==19 && row==142+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==19 && row==162+(ch-23)*2) flag=0;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==19 && row==206+(ch-46)*2) flag=0;
    for(Int_t ch=57; ch<=62; ++ch) if(layer==19 && row==226+(ch-57)*2) flag=0;
    //AGET-1
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==19 && row==119+ch*2)      flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==19 && row==141+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==19 && row==161+(ch-23)*2) flag=1;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==19 && row==205+(ch-46)*2) flag=1;
    for(Int_t ch=57; ch<=63; ++ch) if(layer==19 && row==225+(ch-57)*2) flag=1;
    //AGET-2
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==18 && row==117+ch*2)      flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==18 && row==139+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==18 && row==159+(ch-23)*2) flag=2;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==18 && row==203+(ch-46)*2) flag=2;
    for(Int_t ch=57; ch<=61; ++ch) if(layer==18 && row==223+(ch-57)*2) flag=2;
    //AGET-3
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==18 && row==116+ch*2)      flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==18 && row==138+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==18 && row==158+(ch-23)*2) flag=3;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==18 && row==202+(ch-46)*2) flag=3;
    for(Int_t ch=57; ch<=61; ++ch) if(layer==18 && row==222+(ch-57)*2) flag=3;
    return flag;
  case 20:
    //AGET-0
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==21 && row==ch*2)          flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==21 && row==22+(ch-12)*2)  flag=0;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==21 && row==42+(ch-23)*2)  flag=0;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==21 && row==86+(ch-46)*2)  flag=0;
    for(Int_t ch=57; ch<=61; ++ch) if(layer==21 && row==106+(ch-57)*2) flag=0;
    //AGET-1
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==21 && row==1+ch*2)        flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==21 && row==23+(ch-12)*2)  flag=1;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==21 && row==43+(ch-23)*2)  flag=1;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==21 && row==87+(ch-46)*2)  flag=1;
    for(Int_t ch=57; ch<=61; ++ch) if(layer==21 && row==107+(ch-57)*2) flag=1;
    //AGET-2
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==20 && row==ch*2)          flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==20 && row==22+(ch-12)*2)  flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==20 && row==42+(ch-23)*2)  flag=2;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==20 && row==86+(ch-46)*2)  flag=2;
    for(Int_t ch=57; ch<=64; ++ch) if(layer==20 && row==106+(ch-57)*2) flag=2;
    //AGET-3
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==20 && row==1+ch*2)        flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==20 && row==23+(ch-12)*2)  flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==20 && row==43+(ch-23)*2)  flag=3;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==20 && row==87+(ch-46)*2)  flag=3;
    for(Int_t ch=57; ch<=64; ++ch) if(layer==20 && row==107+(ch-57)*2) flag=3;
    return flag;
  case 21:
    //AGET-0
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==21 && row==117+ch*2)      flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==21 && row==139+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==21 && row==159+(ch-23)*2) flag=0;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==21 && row==203+(ch-46)*2) flag=0;
    for(Int_t ch=57; ch<=61; ++ch) if(layer==21 && row==223+(ch-57)*2) flag=0;
    //AGET-1
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==21 && row==116+ch*2)      flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==21 && row==138+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==21 && row==158+(ch-23)*2) flag=1;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==21 && row==202+(ch-46)*2) flag=1;
    for(Int_t ch=57; ch<=61; ++ch) if(layer==21 && row==222+(ch-57)*2) flag=1;
    //AGET-2
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==20 && row==123+ch*2)      flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==20 && row==145+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==20 && row==165+(ch-23)*2) flag=2;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==20 && row==209+(ch-46)*2) flag=2;
    for(Int_t ch=57; ch<=64; ++ch) if(layer==20 && row==229+(ch-57)*2) flag=2;
    //AGET-3
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==20 && row==122+ch*2)      flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==20 && row==144+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==20 && row==164+(ch-23)*2) flag=3;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==20 && row==208+(ch-46)*2) flag=3;
    for(Int_t ch=57; ch<=64; ++ch) if(layer==20 && row==228+(ch-57)*2) flag=3;
    return flag;
  case 22:
    //AGET-0
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==23 && row==1+ch*2)       flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==23 && row==23+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==23 && row==43+(ch-23)*2) flag=0;
    for(Int_t ch=46; ch<=54; ++ch) if(layer==23 && row==87+(ch-46)*2) flag=0;
    //AGET-1
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==23 && row==ch*2)         flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==23 && row==22+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==23 && row==42+(ch-23)*2) flag=1;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==23 && row==86+(ch-46)*2) flag=1;
    //AGET-2
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==22 && row==1+ch*2)       flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==22 && row==23+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=35; ++ch) if(layer==22 && row==43+(ch-23)*2) flag=2;
    for(Int_t ch=39; ch<=44; ++ch) if(layer==22 && row==75+(ch-39)*2) flag=2;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==22 && row==87+(ch-46)*2) flag=2;
    if(layer==22 && row==107) flag=2;
    //AGET-3
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==22 && row==ch*2)          flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==22 && row==22+(ch-12)*2)  flag=3;
    for(Int_t ch=23; ch<=36; ++ch) if(layer==22 && row==42+(ch-23)*2)  flag=3;
    for(Int_t ch=39; ch<=44; ++ch) if(layer==22 && row==74+(ch-39)*2)  flag=3;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==22 && row==86+(ch-46)*2)  flag=3;
    for(Int_t ch=57; ch<=58; ++ch) if(layer==22 && row==106+(ch-57)*2) flag=3;
    return flag;
  case 23:
    //AGET-0
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==23 && row==106+ch*2)      flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==23 && row==128+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==23 && row==148+(ch-23)*2) flag=0;
    for(Int_t ch=46; ch<=54; ++ch) if(layer==23 && row==192+(ch-46)*2) flag=0;
    //AGET-1
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==23 && row==105+ch*2)      flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==23 && row==127+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==23 && row==147+(ch-23)*2) flag=1;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==23 && row==191+(ch-46)*2) flag=1;
    //AGET-2
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==22 && row==110+ch*2)      flag=2;
    for(Int_t ch=12; ch<=17; ++ch) if(layer==22 && row==132+(ch-12)*2) flag=2;
    if(layer==22 && row==150) flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==22 && row==152+(ch-23)*2) flag=2;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==22 && row==196+(ch-46)*2) flag=2;
    if(layer==22 && row==216) flag=2;
    //AGET-3
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==22 && row==109+ch*2)      flag=3;
    for(Int_t ch=12; ch<=18; ++ch) if(layer==22 && row==131+(ch-12)*2) flag=3;
    if(layer==22 && row==149) flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==22 && row==151+(ch-23)*2) flag=3;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==22 && row==195+(ch-46)*2) flag=3;
    for(Int_t ch=57; ch<=58; ++ch) if(layer==22 && row==215+(ch-57)*2) flag=3;
    return flag;
  case 24:
    //AGET-0
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==25 && row==1+ch*2)       flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==25 && row==23+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==25 && row==43+(ch-23)*2) flag=0;
    for(Int_t ch=46; ch<=52; ++ch) if(layer==25 && row==87+(ch-46)*2) flag=0;
    //AGET-1
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==25 && row==ch*2)         flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==25 && row==22+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==25 && row==42+(ch-23)*2) flag=1;
    for(Int_t ch=46; ch<=53; ++ch) if(layer==25 && row==86+(ch-46)*2) flag=1;
    //AGET-2
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==24 && row==1+ch*2)       flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==24 && row==23+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=28; ++ch) if(layer==24 && row==43+(ch-23)*2) flag=2;
    for(Int_t ch=32; ch<=44; ++ch) if(layer==24 && row==61+(ch-32)*2) flag=2;
    for(Int_t ch=46; ch<=53; ++ch) if(layer==24 && row==87+(ch-46)*2) flag=2;
    //AGET-3
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==24 && row==ch*2)         flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==24 && row==22+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=29; ++ch) if(layer==24 && row==42+(ch-23)*2) flag=3;
    for(Int_t ch=33; ch<=44; ++ch) if(layer==24 && row==62+(ch-33)*2) flag=3;
    for(Int_t ch=46; ch<=54; ++ch) if(layer==24 && row==86+(ch-46)*2) flag=3;
    return flag;
  case 25:
    //AGET-0
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==25 && row==102+ch*2)      flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==25 && row==124+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==25 && row==144+(ch-23)*2) flag=0;
    for(Int_t ch=46; ch<=52; ++ch) if(layer==25 && row==188+(ch-46)*2) flag=0;
    //AGET-1
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==25 && row==101+ch*2)      flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==25 && row==123+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==25 && row==143+(ch-23)*2) flag=1;
    for(Int_t ch=46; ch<=53; ++ch) if(layer==25 && row==187+(ch-46)*2) flag=1;
    //AGET-2
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==24 && row==104+ch*2)      flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==24 && row==126+(ch-12)*2) flag=2;
    for(Int_t ch=26; ch<=44; ++ch) if(layer==24 && row==152+(ch-26)*2) flag=2;
    for(Int_t ch=46; ch<=53; ++ch) if(layer==24 && row==190+(ch-46)*2) flag=2;
    //AGET-3
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==24 && row==103+ch*2)      flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==24 && row==125+(ch-12)*2) flag=3;
    for(Int_t ch=26; ch<=44; ++ch) if(layer==24 && row==151+(ch-26)*2) flag=3;
    for(Int_t ch=46; ch<=54; ++ch) if(layer==24 && row==189+(ch-46)*2) flag=3;
    return flag;
  case 26:
    //AGET-0
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==27 && row==ch*2)         flag=0;
    for(Int_t ch=12; ch<=20; ++ch) if(layer==27 && row==22+(ch-12)*2) flag=0;
    for(Int_t ch=24; ch<=44; ++ch) if(layer==27 && row==44+(ch-24)*2) flag=0;
    for(Int_t ch=46; ch<=51; ++ch) if(layer==27 && row==86+(ch-46)*2) flag=0;
    //AGET-1
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==27 && row==1+ch*2)       flag=1;
    for(Int_t ch=12; ch<=19; ++ch) if(layer==27 && row==23+(ch-12)*2) flag=1;
    if(layer==27 && row==41) flag=1;
    for(Int_t ch=24; ch<=44; ++ch) if(layer==27 && row==45+(ch-24)*2) flag=1;
    for(Int_t ch=46; ch<=51; ++ch) if(layer==27 && row==87+(ch-46)*2) flag=1;
    //AGET-2
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==26 && row==ch*2)         flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==26 && row==22+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==26 && row==42+(ch-23)*2) flag=2;
    for(Int_t ch=46; ch<=52; ++ch) if(layer==26 && row==86+(ch-46)*2) flag=2;
    //AGET-3
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==26 && row==1+ch*2)       flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==26 && row==23+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==26 && row==43+(ch-23)*2) flag=3;
    for(Int_t ch=46; ch<=52; ++ch) if(layer==26 && row==87+(ch-46)*2) flag=3;
    return flag;
  case 27:
    //AGET-0
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==27 && row==99+ch*2)       flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==27 && row==121+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=28; ++ch) if(layer==27 && row==141+(ch-23)*2) flag=0;
    for(Int_t ch=31; ch<=44; ++ch) if(layer==27 && row==157+(ch-31)*2) flag=0;
    for(Int_t ch=46; ch<=51; ++ch) if(layer==27 && row==185+(ch-46)*2) flag=0;
    //AGET-1
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==27 && row==98+ch*2)       flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==27 && row==120+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=28; ++ch) if(layer==27 && row==140+(ch-23)*2) flag=1;
    if(layer==27 && row==154) flag=1;
    for(Int_t ch=32; ch<=44; ++ch) if(layer==27 && row==158+(ch-32)*2) flag=1;
    for(Int_t ch=46; ch<=51; ++ch) if(layer==27 && row==184+(ch-46)*2) flag=1;
    //AGET-2
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==26 && row==101+ch*2)      flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==26 && row==123+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==26 && row==143+(ch-23)*2) flag=2;
    for(Int_t ch=46; ch<=52; ++ch) if(layer==26 && row==187+(ch-46)*2) flag=2;
    //AGET-3
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==26 && row==100+ch*2)      flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==26 && row==122+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==26 && row==142+(ch-23)*2) flag=3;
    for(Int_t ch=46; ch<=52; ++ch) if(layer==26 && row==186+(ch-46)*2) flag=3;
    return flag;
  case 28:
    //AGET-0
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==29 && row==1+ch*2)       flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==29 && row==23+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=33; ++ch) if(layer==29 && row==43+(ch-23)*2) flag=0;
    //AGET-1
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==29 && row==ch*2)         flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==29 && row==22+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=34; ++ch) if(layer==29 && row==42+(ch-23)*2) flag=1;
    //AGET-2
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==28 && row==1+ch*2)       flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==28 && row==23+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==28 && row==43+(ch-23)*2) flag=2;
    if(layer==28 && row==87) flag=2;
    //AGET-3
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==28 && row==ch*2)         flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==28 && row==22+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==28 && row==42+(ch-23)*2) flag=3;
    if(layer==28 && row==86) flag=3;
    if(layer==28 && row==88) flag=3;
    return flag;
  case 29:
    //AGET-0
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==29 && row==66+ch*2)       flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==29 && row==88+(ch-12)*2)  flag=0;
    for(Int_t ch=23; ch<=33; ++ch) if(layer==29 && row==108+(ch-23)*2) flag=0;
    //AGET-1
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==29 && row==65+ch*2)       flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==29 && row==87+(ch-12)*2)  flag=1;
    for(Int_t ch=23; ch<=34; ++ch) if(layer==29 && row==107+(ch-23)*2) flag=1;
    //AGET-2
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==28 && row==90+ch*2)       flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==28 && row==112+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==28 && row==132+(ch-23)*2) flag=2;
    if(layer==28 && row==176) flag=2;
    //AGET-3
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==28 && row==89+ch*2)       flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==28 && row==111+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==28 && row==131+(ch-23)*2) flag=3;
    if(layer==28 && row==175) flag=3;
    if(layer==28 && row==177) flag=3;
    return flag;
  case 30:
    //AGET-0
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==31 && row==1+ch*2)       flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==31 && row==23+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==31 && row==43+(ch-23)*2) flag=0;
    if(layer==31 && row==87) flag=0;
    if(layer==31 && row==89) flag=0;
    //AGET-1
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==31 && row==ch*2)         flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==31 && row==22+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==31 && row==42+(ch-23)*2) flag=1;
    if(layer==31 && row==86) flag=1;
    if(layer==31 && row==88) flag=1;
    //AGET-2
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==30 && row==ch*2)         flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==30 && row==22+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==30 && row==42+(ch-23)*2) flag=2;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==30 && row==86+(ch-46)*2) flag=2;
    if(layer==30 && row==106) flag=2;
    //AGET-3
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==30 && row==1+ch*2)       flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==30 && row==23+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==30 && row==43+(ch-23)*2) flag=3;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==30 && row==87+(ch-46)*2) flag=3;
    if(layer==30 && row==107) flag=3;
    return flag;
  default:
    return -1;
  }
}

//_____________________________________________________________________________
TpcPadParam*
TpcPadHelper::GetParam( Int_t asad, Int_t aget, Int_t channel ) const
{
  Int_t key = TpcPadParam::Key( asad, aget, channel );
  TpcPadParam* param = nullptr;
  auto itr = m_map.find( key );
  if( itr != m_map.end() )
    param = itr->second;
#if 0
  else
    std::cerr << FUNC_NAME << " No parameter, AsAd=" << asad
	      << ", AGet=" << aget << ", Channel=" << channel << std::endl;
#endif
  return param;
}

//_____________________________________________________________________________
TpcPadParam*
TpcPadHelper::GetParam( Int_t layer, Int_t row ) const
{
  for( const auto& p : m_map ){
    if( p.second->LayerId() == layer &&
        p.second->RowId() == row ){
      return p.second;
    }
  }
#if 0
  std::cerr << FUNC_NAME << " No parameter, Layer=" << layer
            << ", Row=" << row << std::endl;
#endif
  return nullptr;
}

//_____________________________________________________________________________
TpcPadParam*
TpcPadHelper::GetParam( Int_t pad ) const
{
  for( const auto& p : m_map ){
    if( p.second->PadId() == pad )
      return p.second;
  }
#if 0
  std::cerr << FUNC_NAME << " No parameter, Pad=" << pad << std::endl;
#endif
  return nullptr;
}

//_____________________________________________________________________________
TVector3
TpcPadHelper::GetPoint( Int_t pad ) const
{
  TVector3 vec;
  Int_t layer, row;
  Int_t sum=0;
  for( layer=0; layer<=30 && sum+TpcPadHelper::PadParameter[layer][1]<=pad; layer++ ){
    sum += TpcPadHelper::PadParameter[layer][1];
  }
  row = pad - sum;

  if( row>TpcPadHelper::PadParameter[layer][1] ) {
    vec.SetXYZ(0,-1,0);
  } else {
    Double_t stheta = 180.-(360./TpcPadHelper::PadParameter[layer][3])*TpcPadHelper::PadParameter[layer][1]/2.;
    Double_t theta = stheta+(row+0.5)*360./TpcPadHelper::PadParameter[layer][3]-180.;

    Double_t x = TpcPadHelper::PadParameter[layer][2]*sin(theta*acos(-1)/180.);
    Double_t z = TpcPadHelper::PadParameter[layer][2]*cos(theta*acos(-1)/180.)-143.;
    vec.SetXYZ(x,0,z);
  }
  return vec;
}

//______________________________________________________________________________
void
ConfMan::initializeTpcPadHelper()
{
  if(name_file_["TPCPAD:"] != ""){
    TpcPadHelper& gTpcPad = TpcPadHelper::GetInstance();
    gTpcPad.SetFileName(name_file_["TPCPAD:"]);
    flag_[kIsGood] = gTpcPad.Initialize();
  }else{
    std::cout << "#E ConfMan::"
	      << " File path does not exist in " << name_file_["TPCPAD:"] 
	      << std::endl;
    flag_.reset(kIsGood);
  }
}

//_____________________________________________________________________________
/*
Bool_t
TpcPadHelper::Initialize( const TString &file_name )
{
  std::ifstream ifs( file_name );
  if( !ifs.is_open() ){
    std::cerr << FUNC_NAME << " " << file_name
	      << " ... failed to open" << std::endl;
    return false;
  }
  TString line;
  while( ifs.good() && line.ReadLine( ifs ) ){
    if( line[0]=='#' )
      continue;
    std::istringstream iss( line.Data() );
    Int_t aget, asad, channel, layer, row, pad;
    iss >> aget >> asad >> channel >> layer >> row >> pad;
    if( aget < 0 || asad < 0 || channel < 0 )
      continue;
    Int_t key = TpcPadParam::Key( asad, aget, channel );
    // std::cout << line << std::endl;
    if( m_map.find( key ) != m_map.end() ){
      std::cerr << FUNC_NAME << " " << file_name
		<< " ... duplicated parameter" << std::endl;
      TpcPadParam( asad, aget, channel, pad, layer, row ).Print();
      m_map.find( key )->second->Print();
      continue;
      // return false;
    }
    if( GetParam( pad ) ){
      std::cerr << FUNC_NAME << " " << file_name
		<< " ... duplicated parameter" << std::endl;
      TpcPadParam( asad, aget, channel, pad, layer, row ).Print();
      GetParam( pad )->Print();
      continue;
      // return false;
    }
    m_map[key] = new TpcPadParam( asad, aget, channel,
                                  pad, layer, row );
  }
#if 0
  std::cout << FUNC_NAME << " " << file_name
	    << " ... initialized (NParam=" << m_map.size() << ")" << std::endl;
  for( const auto& p : m_map ){
    std::cout << " PadId=" << p.second->PadId()
	      << " LayerId=" << p.second->LayerId()
	      << " RowId=" << p.second->RowId()
	      << std::endl;
  }
#endif
  return true;
}
*/


Bool_t
TpcPadHelper::Initialize( const std::string& filename )
{
  m_file_name = filename;
  return Initialize();
};

Bool_t
TpcPadHelper::Initialize( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  std::ifstream ifs( m_file_name.c_str() );
  if( !ifs.is_open() ){
    std::cerr << "#E " << func_name << " "
		<< "No such parameter file : " << m_file_name << std::endl;
    std::exit(EXIT_FAILURE);
  }
  
  TString line;
  while( ifs.good() && line.ReadLine( ifs ) ){
    if( line[0]=='#' )
      continue;
    std::istringstream iss( line.Data() );
    Int_t aget, asad, channel, layer, row, pad;
    iss >> aget >> asad >> channel >> layer >> row >> pad;
    if( aget < 0 || asad < 0 || channel < 0 )
      continue;
    Int_t key = TpcPadParam::Key( asad, aget, channel );
    if( m_map.find( key ) != m_map.end() ){
      std::cerr << FUNC_NAME << " " << m_file_name
		<< " ... duplicated parameter" << std::endl;
      TpcPadParam( asad, aget, channel, pad, layer, row ).Print();
      m_map.find( key )->second->Print();
      continue;
      // return false;
    }
    if( GetParam( pad ) ){
      std::cerr << FUNC_NAME << " " << m_file_name
		<< " ... duplicated parameter" << std::endl;
      TpcPadParam( asad, aget, channel, pad, layer, row ).Print();
      GetParam( pad )->Print();
      continue;
      // return false;
    }
    m_map[key] = new TpcPadParam( asad, aget, channel,
                                  pad, layer, row );
  }
#if 0
  std::cout << FUNC_NAME << " " << m_file_name
	    << " ... initialized (NParam=" << m_map.size() << ")" << std::endl;
  for( const auto& p : m_map ){
    std::cout << " PadId=" << p.second->PadId()
	      << " LayerId=" << p.second->LayerId()
	      << " RowId=" << p.second->RowId()
	      << std::endl;
  }
#endif
  return true;
}
