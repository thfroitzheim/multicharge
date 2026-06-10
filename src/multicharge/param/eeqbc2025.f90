! This file is part of multicharge.
! SPDX-Identifier: Apache-2.0
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

!> Bond capacitor electronegativity equilibration charge model published in
!>
!> Thomas Froitzheim, Marcel Müller, Andreas Hansen, and Stefan Grimme,
!> *J. Chem. Phys.*, **2025**, 162, 214109.
!> DOI: [10.1063/5.0268978](https://dx.doi.org/10.1063/5.0268978)
!> 
!> Updated from of the parametrization and minor model changes published in
!>
!> Thomas Froitzheim, Marcel Müller, Andreas Hansen, and Stefan Grimme,
!> *ChemRxiv*, **2025**.
!> DOI: [10.26434/chemrxiv-2025-bjxvt](https://doi.org/10.26434/chemrxiv-2025-bjxvt)
!>
!> The original parametrization can be used in multicharge v0.5.0. 
module multicharge_param_eeqbc2025
   use mctc_env, only: wp
   use mctc_io_symbols, only: to_number
   implicit none
   private

   public :: get_eeqbc_chi, get_eeqbc_eta, get_eeqbc_rad, get_eeqbc_kcnchi, &
      & get_eeqbc_kqchi, get_eeqbc_kqeta, get_eeqbc_cap, &
      & get_eeqbc_cov_radii, get_eeqbc_avg_cn, get_eeqbc_rvdw_scale

   !> Element-specific electronegativity for the EEQ_BC charges.
   interface get_eeqbc_chi
      module procedure :: get_eeqbc_chi_sym
      module procedure :: get_eeqbc_chi_num
   end interface get_eeqbc_chi

   !> Element-specific chemical hardnesses for the EEQ_BC charges.
   interface get_eeqbc_eta
      module procedure :: get_eeqbc_eta_sym
      module procedure :: get_eeqbc_eta_num
   end interface get_eeqbc_eta

   !> Element-specific charge widths for the EEQ_BC charges.
   interface get_eeqbc_rad
      module procedure :: get_eeqbc_rad_sym
      module procedure :: get_eeqbc_rad_num
   end interface get_eeqbc_rad

   !> Element-specific CN scaling of the electronegativity for the EEQ_BC charges.
   interface get_eeqbc_kcnchi
      module procedure :: get_eeqbc_kcnchi_sym
      module procedure :: get_eeqbc_kcnchi_num
   end interface get_eeqbc_kcnchi

   !> Element-specific local q scaling of the electronegativity for the EEQ_BC charges.
   interface get_eeqbc_kqchi
      module procedure :: get_eeqbc_kqchi_sym
      module procedure :: get_eeqbc_kqchi_num
   end interface get_eeqbc_kqchi

   !> Element-specific local q scaling of the chemical hardness for the EEQ_BC charges.
   interface get_eeqbc_kqeta
      module procedure :: get_eeqbc_kqeta_sym
      module procedure :: get_eeqbc_kqeta_num
   end interface get_eeqbc_kqeta

   !> Element-specific bond capacitance for the EEQ_BC charges.
   interface get_eeqbc_cap
      module procedure :: get_eeqbc_cap_sym
      module procedure :: get_eeqbc_cap_num
   end interface get_eeqbc_cap

   !> Element-specific covalent radii for the CN for the EEQ_BC charges.
   interface get_eeqbc_cov_radii
      module procedure :: get_eeqbc_cov_radii_sym
      module procedure :: get_eeqbc_cov_radii_num
   end interface get_eeqbc_cov_radii

   !> Element-specific average CN for the EEQ_BC charges.
   interface get_eeqbc_avg_cn
      module procedure :: get_eeqbc_avg_cn_sym
      module procedure :: get_eeqbc_avg_cn_num
   end interface get_eeqbc_avg_cn

   !> Element-specific scaling of the pairwise van-der-Waals radii.
   interface get_eeqbc_rvdw_scale
      module procedure :: get_eeqbc_rvdw_scale_sym
      module procedure :: get_eeqbc_rvdw_scale_num
   end interface get_eeqbc_rvdw_scale


   !> Maximum atomic number allowed in EEQ_BC calculations
   integer, parameter :: max_elem = 103

   !> Element-specific electronegativity for the EEQ_BC charges.
real(wp), parameter :: eeqbc_chi(max_elem) = [&
      &  0.9586845058_wp, -0.1832585081_wp, -0.1994511709_wp,  0.6274825494_wp, & !1-4
      &  0.9097045017_wp,  1.8344318686_wp,  2.6850970494_wp,  3.0506357235_wp, & !5-8
      &  2.5120529912_wp,  0.7163540834_wp, -1.0908240913_wp, -0.0779944731_wp, & !9-12
      &  0.0522075411_wp,  0.8262209147_wp,  1.4502301986_wp,  2.2074275540_wp, & !13-16
      &  2.5018284543_wp,  0.3699716844_wp, -2.0130931696_wp, -0.3192109028_wp, & !17-20
      & -0.5751500114_wp, -0.0605203088_wp,  0.0404113075_wp,  0.1541675079_wp, & !21-24
      &  0.0959020251_wp,  0.1629816983_wp,  0.0721490965_wp,  0.5995588046_wp, & !25-28
      &  0.2889279545_wp, -0.0155186655_wp,  0.1705506933_wp,  0.8837385035_wp, & !29-32
      &  1.4252146278_wp,  2.6998652654_wp,  1.8959112567_wp,  0.2704832340_wp, & !33-36
      & -1.2978905925_wp, -0.6482447696_wp, -0.3989214011_wp, -0.1737399062_wp, & !37-40
      & -0.1744722135_wp,  0.4200884618_wp,  0.5517684435_wp,  0.3448618275_wp, & !41-44
      &  0.3181241233_wp,  0.9943064927_wp, -0.0280463867_wp, -0.0848711274_wp, & !45-48
      &  0.1103620247_wp,  0.7360044259_wp,  0.8092608583_wp,  1.7802168727_wp, & !49-52
      &  1.7270720040_wp,  0.2260612052_wp, -1.5680486332_wp, -0.6590352941_wp, & !53-56
      & -0.9060549690_wp, -1.3316329034_wp, -1.8150039750_wp, -1.5011820425_wp, & !57-60
      & -1.7412826677_wp, -2.0563736869_wp, -1.7649191046_wp, -1.5001192705_wp, & !61-64
      & -1.6096777876_wp, -1.5270862651_wp, -1.2949666810_wp, -1.7415737194_wp, & !65-68
      & -1.6761373483_wp, -1.5483210565_wp, -1.4223252423_wp, -0.2374493733_wp, & !69-72
      &  0.2801493770_wp,  0.5767462052_wp,  0.9471893926_wp,  0.8361727489_wp, & !73-76
      &  0.7144855704_wp,  1.1569981859_wp,  0.4635644309_wp,  0.2201241678_wp, & !77-80
      &  0.2470285233_wp,  0.0937218923_wp,  0.5204647797_wp,  1.3096619310_wp, & !81-84
      &  1.5362656547_wp,  0.0723900534_wp, -1.5101245871_wp, -0.4456908449_wp, & !85-88
      & -0.7883619246_wp, -0.6215317730_wp, -0.8978780690_wp, -1.0218594922_wp, & !89-92
      & -0.5669120729_wp, -0.7828450851_wp, -0.7240270747_wp, -0.7511566794_wp, & !93-96
      & -0.5594997284_wp, -0.5617146680_wp, -0.7604827875_wp, -0.9715813282_wp, & !97-100
      & -0.7720799645_wp, -0.9248211063_wp, -1.0394830497_wp] !101-103

   !> Element-specific chemical hardnesses for the EEQ_BC charges.
   real(wp), parameter :: eeqbc_eta(max_elem) = [&
      &  6.2117434795_wp,  1.0037450542_wp, -0.3072745674_wp, -0.1655235755_wp, & !1-4
      & -0.7065327502_wp, -0.8905032667_wp, -0.9667138384_wp, -0.5302539432_wp, & !5-8
      & -2.8921939301_wp, -0.8713602610_wp, -2.7449394352_wp, -2.4721871938_wp, & !9-12
      & -2.5822029727_wp, -1.5706776705_wp, -1.7175610731_wp, -0.9211245787_wp, & !13-16
      & -0.7311938938_wp, -2.2082005781_wp,  0.8539343179_wp, -0.4736934466_wp, & !17-20
      & -0.4668674639_wp, -0.5635709055_wp, -0.9477885754_wp, -0.4724139899_wp, & !21-24
      & -0.7802954666_wp, -1.3114434237_wp, -1.2804825781_wp, -0.3095352927_wp, & !25-28
      & -0.2239635680_wp, -0.3180242538_wp, -0.3412072442_wp, -1.4713683540_wp, & !29-32
      & -1.3048043637_wp, -0.3445602696_wp,  0.9657440120_wp,  1.0614087719_wp, & !33-36
      &  0.5135338808_wp, -1.9743278013_wp, -1.2578555157_wp, -0.6463997763_wp, & !37-40
      & -0.1343396385_wp,  0.0817323135_wp,  0.1405555844_wp,  0.2283011469_wp, & !41-44
      &  0.2622937353_wp,  0.0972819315_wp,  0.4736046942_wp, -1.5228469757_wp, & !45-48
      & -1.2235567223_wp, -1.6338397262_wp, -1.5316779125_wp, -0.3361857147_wp, & !49-52
      &  0.9717354739_wp,  1.0610988738_wp,  0.4556834722_wp, -1.1771590144_wp, & !53-56
      & -1.2305197024_wp, -2.2734199694_wp, -1.4809108807_wp, -1.5112401585_wp, & !57-60
      & -1.5903079899_wp, -1.8149463240_wp, -1.7060358119_wp, -1.4059798146_wp, & !61-64
      & -1.9733375334_wp, -2.4546649588_wp, -1.9064440299_wp, -1.6661613206_wp, & !65-68
      & -0.5606406060_wp, -0.6628227474_wp, -0.5906274026_wp, -0.6825996599_wp, & !69-72
      & -0.4600687602_wp, -0.1457719211_wp,  0.1495310516_wp,  0.3395136550_wp, & !73-76
      &  0.2783910880_wp,  0.2294794780_wp,  0.1979502960_wp,  0.1019963149_wp, & !77-80
      & -0.2387064009_wp, -0.4138016801_wp, -1.2738664162_wp, -0.9936574116_wp, & !81-84
      & -0.7277636729_wp,  0.2585497078_wp, -0.1286070329_wp, -3.4669302825_wp, & !85-88
      & -0.8619111416_wp, -0.8481567043_wp, -0.7426616457_wp, -0.5841985797_wp, & !89-92
      & -0.5503077938_wp, -0.5350180167_wp, -0.5270374907_wp, -0.4286554958_wp, & !93-96
      & -0.6115218147_wp, -1.0257953048_wp, -0.9333857527_wp, -0.7285519539_wp, & !97-100
      & -1.0341236820_wp, -1.5525834066_wp, -2.3048144101_wp] !101-103

   !> Element-specific charge widths for the EEQ_BC charges.
   real(wp), parameter :: eeqbc_rad(max_elem) = [&
      &  5.5912707296_wp,  0.0319680252_wp,  0.2654744914_wp,  0.4783453222_wp, & !1-4
      &  0.3136998131_wp,  0.2662238318_wp,  0.2345269832_wp,  0.2831846847_wp, & !5-8
      &  0.1578455294_wp,  0.0379636282_wp,  0.1531318295_wp,  0.2144817128_wp, & !9-12
      &  0.1751357840_wp,  0.2498032387_wp,  0.2260990532_wp,  0.2590540342_wp, & !13-16
      &  0.2368013898_wp,  0.1553109973_wp,  0.6118652728_wp,  0.4931516154_wp, & !17-20
      &  0.4197132100_wp,  0.4273555097_wp,  0.3365100687_wp,  0.4483070645_wp, & !21-24
      &  0.4644391977_wp,  0.3695472354_wp,  0.3023717098_wp,  0.5290574593_wp, & !25-28
      &  0.4048963513_wp,  0.4140755332_wp,  0.3345108417_wp,  0.2850759440_wp, & !29-32
      &  0.3023889721_wp,  0.2920352851_wp,  0.5709767735_wp,  0.6403833876_wp, & !33-36
      &  0.4617736621_wp,  0.2791647245_wp,  0.3642855711_wp,  0.4081368797_wp, & !37-40
      &  0.4503724504_wp,  0.6392739006_wp,  0.8352943116_wp,  0.5833279825_wp, & !41-44
      &  0.6963325865_wp,  0.7371968337_wp,  0.4324434219_wp,  0.2620697372_wp, & !45-48
      &  0.3099888904_wp,  0.2821904756_wp,  0.2887594983_wp,  0.4266667003_wp, & !49-52
      &  0.6363726306_wp,  0.7944750000_wp,  0.5395040208_wp,  0.4267973823_wp, & !53-56
      &  0.3798639638_wp,  0.2260316791_wp,  0.2639463325_wp,  0.2523397692_wp, & !57-60
      &  0.2572790476_wp,  0.2495919396_wp,  0.2370655997_wp,  0.2730743122_wp, & !61-64
      &  0.2304725966_wp,  0.2073484814_wp,  0.2255115944_wp,  0.2347722551_wp, & !65-68
      &  0.3144718549_wp,  0.2837958535_wp,  0.3317696869_wp,  0.4209722417_wp, & !69-72
      &  0.5258011341_wp,  0.6950629571_wp,  0.9441290032_wp,  0.9350134020_wp, & !73-76
      &  0.8423443858_wp,  0.6296021058_wp,  0.5352064650_wp,  0.4754404541_wp, & !77-80
      &  0.4249487814_wp,  0.5071710405_wp,  0.3408674091_wp,  0.3543727341_wp, & !81-84
      &  0.3972360834_wp,  0.4947284375_wp,  0.2245570124_wp,  0.1860151688_wp, & !85-88
      &  0.4707814851_wp,  0.5201838290_wp,  0.4953693956_wp,  0.5141120022_wp, & !89-92
      &  0.5816032792_wp,  0.5073515687_wp,  0.5087851219_wp,  0.5265267914_wp, & !93-96
      &  0.4788745924_wp,  0.3704594581_wp,  0.3434101217_wp,  0.4258264048_wp, & !97-100
      &  0.3339586119_wp,  0.2588437838_wp,  0.2231145137_wp] !101-103

   !> Element-specific CN scaling of the electronegativity for the EEQ_BC charges.
   real(wp), parameter :: eeqbc_kcnchi(max_elem) = [&
      &  6.9030913107_wp, 10.8905314394_wp, -1.7295644750_wp, -0.3499322614_wp, & !1-4
      &  0.3858214020_wp,  1.3491730337_wp,  1.7329254379_wp,  2.7058472481_wp, & !5-8
      &  3.9622085001_wp, 15.2706951951_wp, -1.9327051470_wp, -0.8693519203_wp, & !9-12
      &  3.0762641816_wp,  0.9067547703_wp,  1.1927952927_wp,  2.0888194707_wp, & !13-16
      &  5.8116044898_wp, 23.1012659552_wp, -4.7731233351_wp, -2.3748625547_wp, & !17-20
      & -1.5289932630_wp, -1.8950611938_wp, -1.6065357296_wp, -2.2311613272_wp, & !21-24
      & -2.4395323633_wp, -2.0914311293_wp, -1.7379620595_wp, -1.4850174452_wp, & !25-28
      & -0.4528169887_wp, -0.4144566585_wp,  0.0904944221_wp,  0.0575852483_wp, & !29-32
      &  0.9851764195_wp,  3.5480797118_wp,  3.5825743822_wp,  6.3334974446_wp, & !33-36
      & -2.1608751327_wp, -2.0170138692_wp, -1.4025251857_wp, -1.5333855409_wp, & !37-40
      & -1.2010276307_wp, -0.5452897086_wp, -1.3484357806_wp, -1.0835634306_wp, & !41-44
      & -0.9183133691_wp, -0.0510577333_wp,  0.1668945162_wp,  0.6870176069_wp, & !45-48
      &  0.2978404529_wp,  0.6423075130_wp,  0.5305891471_wp,  2.1182567221_wp, & !49-52
      &  3.7381642025_wp,  5.6080950757_wp, -6.1564029272_wp, -3.2554261282_wp, & !53-56
      & -2.2540025014_wp, -1.3083435188_wp, -2.4408676628_wp, -2.1847268730_wp, & !57-60
      & -2.6273861457_wp, -2.4951121649_wp, -1.9994928747_wp, -2.1904532197_wp, & !61-64
      & -2.2249068931_wp, -2.2018620294_wp, -1.9251039162_wp, -1.8091572806_wp, & !65-68
      & -1.5622167680_wp, -1.5714049267_wp, -1.3909397007_wp, -2.3204453180_wp, & !69-72
      & -1.2444422792_wp, -0.3112764249_wp, -0.6803075093_wp, -0.7778572521_wp, & !73-76
      & -0.7271526668_wp, -0.0694782588_wp,  0.8489527334_wp,  0.0801883948_wp, & !77-80
      &  0.0733651637_wp,  0.0812227597_wp,  0.0826816800_wp,  1.5286499679_wp, & !81-84
      &  2.7557266288_wp,  6.2345982648_wp, -6.0197434166_wp, -2.9008271371_wp, & !85-88
      & -5.7117540178_wp, -8.1840504899_wp, -4.2745612221_wp, -3.7495812060_wp, & !89-92
      & -3.5528875479_wp, -2.6067625276_wp, -3.1912184276_wp, -2.8781517067_wp, & !93-96
      & -1.8144214047_wp, -1.8529475361_wp, -1.9364462087_wp, -3.1496116814_wp, & !97-100
      & -1.8455752387_wp, -1.1047153940_wp, -1.1443332524_wp] !101-103

   !> Element-specific local q scaling of the electronegativity for the EEQ_BC charges.
   real(wp), parameter :: eeqbc_kqchi(max_elem) = [&
      &  2.6206470553_wp,  0.1770760012_wp, 10.5972942390_wp,  7.5027910577_wp, & !1-4
      &  7.1446422226_wp,  4.8217084331_wp,  5.2003647626_wp,  4.7549734519_wp, & !5-8
      &  4.9568179848_wp,  0.3214572331_wp, 11.5147437076_wp,  8.8223615912_wp, & !9-12
      &  8.3224185324_wp,  8.4903380327_wp,  6.9724776594_wp,  6.8674536810_wp, & !13-16
      &  6.4316079018_wp,  0.8524422257_wp, 10.9309517877_wp,  7.9066882566_wp, & !17-20
      &  7.9082964972_wp,  8.1123014845_wp,  6.0302357804_wp,  6.5614825613_wp, & !21-24
      &  5.5064263601_wp,  6.3377764464_wp,  7.4024378919_wp,  7.5937830592_wp, & !25-28
      &  6.7107692580_wp,  7.4459004777_wp,  8.9706608738_wp,  8.1462029165_wp, & !29-32
      &  7.5448879697_wp,  7.9410516980_wp,  7.3325921931_wp,  1.9349444137_wp, & !33-36
      &  9.6468730632_wp,  8.4824325317_wp,  7.9660136942_wp,  8.1552227366_wp, & !37-40
      &  7.2782011535_wp,  5.9653455492_wp,  6.1895410568_wp,  6.8888595191_wp, & !41-44
      &  6.4224580378_wp,  6.7247593062_wp,  6.8703488661_wp,  7.5098829430_wp, & !45-48
      &  8.3783839713_wp,  8.4391966599_wp,  8.6719297344_wp,  7.2080695899_wp, & !49-52
      &  8.2741073157_wp,  2.5294186716_wp, 11.4345707576_wp,  9.1220078938_wp, & !53-56
      &  7.8268762760_wp,  7.5018006411_wp,  6.5487026831_wp,  7.1911843044_wp, & !57-60
      &  7.7013606960_wp,  6.4466301599_wp,  7.1845237593_wp,  7.3587581672_wp, & !61-64
      &  6.8186168534_wp,  7.7363711897_wp,  8.2570260636_wp,  7.3972745267_wp, & !65-68
      &  8.3003305023_wp,  9.0236267802_wp,  8.4346277171_wp,  7.6547683619_wp, & !69-72
      &  7.6868220685_wp,  6.6259959208_wp,  5.7473829947_wp,  6.5622285847_wp, & !73-76
      &  6.7982497723_wp,  7.4464623864_wp,  7.4362192799_wp,  7.2168098334_wp, & !77-80
      &  8.8997589184_wp,  8.6340882085_wp,  7.8544055451_wp,  8.3469191596_wp, & !81-84
      &  7.4739909298_wp,  2.5403077352_wp,  9.0053335516_wp, 10.2266188542_wp, & !85-88
      &  6.9165673940_wp,  6.1504388765_wp,  6.3828633037_wp,  6.0251668716_wp, & !89-92
      &  7.5981451390_wp,  7.5638322247_wp,  7.5272042243_wp,  6.9940864603_wp, & !93-96
      &  7.3532887716_wp,  8.7637119219_wp,  8.1413359913_wp,  7.7539418662_wp, & !97-100
      &  8.0264973298_wp,  7.7905100162_wp,  7.4124232225_wp] !101-103

   !> Element-specific local q scaling of the chemical hardness for the EEQ_BC charges
   real(wp), parameter :: eeqbc_kqeta(max_elem) = [&
      &  7.1684067887_wp, -2.1330695329_wp,  6.1461685618_wp,  2.7902848038_wp, & !1-4
      &  1.8800889740_wp,  0.4935613538_wp,  0.3862385113_wp,  0.5136813591_wp, & !5-8
      &  0.4759394321_wp, -4.2419696886_wp,  9.4972082802_wp,  2.0557995004_wp, & !9-12
      &  1.9845871990_wp,  1.3346327251_wp,  1.2908526903_wp,  0.3416670576_wp, & !13-16
      &  0.3180828436_wp, -3.2993097578_wp, 11.3415213751_wp,  0.7310939211_wp, & !17-20
      &  0.6211282844_wp,  1.2759886000_wp,  0.1921449329_wp,  0.0869508261_wp, & !21-24
      &  0.1515631882_wp,  0.7546094843_wp,  0.8080025000_wp,  0.7285505813_wp, & !25-28
      &  0.9130291798_wp,  0.9645667823_wp,  0.9405434610_wp,  0.1072389701_wp, & !29-32
      &  1.9595103742_wp,  1.1524049835_wp,  1.4620669518_wp, -2.5852384371_wp, & !33-36
      &  4.8449887912_wp,  1.5358189755_wp,  0.5951123126_wp,  0.7249897162_wp, & !37-40
      &  0.1700735025_wp,  0.3015909153_wp,  0.7755351991_wp,  0.7193447842_wp, & !41-44
      &  0.4488641587_wp,  0.2388026807_wp,  0.3801660308_wp,  1.0318081602_wp, & !45-48
      &  1.2114072951_wp,  1.0996763129_wp,  2.9257965601_wp,  0.8099896016_wp, & !49-52
      &  2.4182004117_wp, -2.8233176467_wp,  7.5008417421_wp,  1.5729090621_wp, & !53-56
      &  0.4256819131_wp,  0.7883186408_wp,  0.0288160250_wp,  0.2527020883_wp, & !57-60
      &  0.4032374809_wp,  0.5103582851_wp,  0.4033624033_wp,  0.3760990822_wp, & !61-64
      &  0.3280788122_wp,  0.3668122005_wp,  0.3835627237_wp,  0.4229963913_wp, & !65-68
      &  0.9565235975_wp,  0.9565235975_wp,  0.9565235975_wp,  0.1644465781_wp, & !69-72
      &  0.7534098478_wp,  0.7540761833_wp,  0.3348143039_wp,  0.3235412373_wp, & !73-76
      &  1.0126534620_wp,  0.9302428783_wp,  0.8863455666_wp,  1.4769395922_wp, & !77-80
      &  1.1369484599_wp,  1.7605526375_wp,  1.6774638879_wp,  0.9782289982_wp, & !81-84
      &  1.1037362718_wp, -2.6371025326_wp, -4.2925524576_wp,  4.5876137151_wp, & !85-88
      & -0.1623063675_wp, -0.2973424600_wp, -0.1890013590_wp, -0.4036103842_wp, & !89-92
      &  0.1709685478_wp,  0.8792501286_wp,  0.5978700207_wp,  0.1078659343_wp, & !93-96
      &  0.2052653851_wp,  1.8853241472_wp,  0.6610501164_wp,  0.3532101308_wp, & !97-100
      &  0.5718558621_wp,  0.3646029972_wp,  0.1064889911_wp] !101-103
! &  3.5842033943_wp, -1.0665347665_wp,  3.0730842809_wp,  1.3951424019_wp, & !1-4
! &  0.9400444870_wp,  0.2467806769_wp,  0.1931192556_wp,  0.2568406796_wp, & !5-8
! &  0.2379697161_wp, -2.1209848443_wp,  4.7486041401_wp,  1.0278997502_wp, & !9-12
! &  0.9922935995_wp,  0.6673163626_wp,  0.6454263451_wp,  0.1708335288_wp, & !13-16
! &  0.1590414218_wp, -1.6496548789_wp,  5.6707606875_wp,  0.3655469606_wp, & !17-20
! &  0.3105641422_wp,  0.6379943000_wp,  0.0960724665_wp,  0.0434754131_wp, & !21-24
! &  0.0757815941_wp,  0.3773047421_wp,  0.4040012500_wp,  0.3642752907_wp, & !25-28
! &  0.4565145899_wp,  0.4822833912_wp,  0.4702717305_wp,  0.0536194850_wp, & !29-32
! &  0.9797551871_wp,  0.5762024918_wp,  0.7310334759_wp, -1.2926192186_wp, & !33-36
! &  2.4224943956_wp,  0.7679094878_wp,  0.2975561563_wp,  0.3624948581_wp, & !37-40
! &  0.0850367513_wp,  0.1507954576_wp,  0.3877675996_wp,  0.3596723921_wp, & !41-44
! &  0.2244320794_wp,  0.1194013403_wp,  0.1900830154_wp,  0.5159040801_wp, & !45-48
! &  0.6057036475_wp,  0.5498381565_wp,  1.4628982801_wp,  0.4049948008_wp, & !49-52
! &  1.2091002058_wp, -1.4116588234_wp,  3.7504208711_wp,  0.7864545310_wp, & !53-56
! &  0.2128409565_wp,  0.3941593204_wp,  0.0144080125_wp,  0.1263510442_wp, & !57-60
! &  0.2016187405_wp,  0.2551791426_wp,  0.2016812017_wp,  0.1880495411_wp, & !61-64
! &  0.1640394061_wp,  0.1834061003_wp,  0.1917813618_wp,  0.2114981956_wp, & !65-68
! &  0.4782617988_wp,  0.4782617988_wp,  0.4782617988_wp,  0.0822232890_wp, & !69-72
! &  0.3767049239_wp,  0.3770380917_wp,  0.1674071519_wp,  0.1617706187_wp, & !73-76
! &  0.5063267310_wp,  0.4651214391_wp,  0.4431727833_wp,  0.7384697961_wp, & !77-80
! &  0.5684742299_wp,  0.8802763188_wp,  0.8387319439_wp,  0.4891144991_wp, & !81-84
! &  0.5518681359_wp, -1.3185512663_wp, -2.1462762288_wp,  2.2938068575_wp, & !85-88
! & -0.0811531837_wp, -0.1486712300_wp, -0.0945006795_wp, -0.2018051921_wp, & !89-92
! &  0.0854842739_wp,  0.4396250643_wp,  0.2989350103_wp,  0.0539329671_wp, & !93-96
! &  0.1026326926_wp,  0.9426620736_wp,  0.3305250582_wp,  0.1766050654_wp, & !97-100
! &  0.2859279310_wp,  0.1823014986_wp,  0.0532444955_wp] !101-103

   !> Element-specific bond capacitance for the EEQ_BC charges.
   real(wp), parameter :: eeqbc_cap(max_elem) = [&
      &  0.2538804328_wp,  0.1333013493_wp,  0.3297035680_wp,  0.7624886399_wp, & !1-4
      &  0.6724009815_wp,  1.3833781469_wp,  0.8806944388_wp,  0.8939791850_wp, & !5-8
      &  1.3812298095_wp,  0.5923957985_wp,  0.2618401905_wp,  0.4588107295_wp, & !9-12
      &  0.2392228783_wp,  0.4438090663_wp,  0.6569843019_wp,  0.6518285818_wp, & !13-16
      &  0.8340140753_wp,  0.2437190068_wp,  0.2825077861_wp,  0.8919998262_wp, & !17-20
      &  0.5470331665_wp,  0.6106228196_wp,  0.6255974955_wp,  1.4890183199_wp, & !21-24
      &  1.1408422486_wp,  1.3171844141_wp,  1.2271080175_wp,  1.0092218806_wp, & !25-28
      &  0.9406797225_wp,  0.5799608246_wp,  0.4582522958_wp,  0.5358246591_wp, & !29-32
      &  0.8817039526_wp,  0.5974538308_wp,  1.6537545362_wp,  0.2898446408_wp, & !33-36
      &  0.6389546722_wp,  0.6754089605_wp,  0.7700672811_wp,  0.3975173525_wp, & !37-40
      &  0.6548393363_wp,  1.2407265056_wp,  1.3061976078_wp,  1.3841285711_wp, & !41-44
      &  1.6397584621_wp,  2.1229476499_wp,  0.7005497839_wp,  0.4813303355_wp, & !45-48
      &  0.6142008073_wp,  0.7671354121_wp,  0.7785946357_wp,  0.7109508841_wp, & !49-52
      &  1.7334321442_wp,  0.3312362125_wp,  0.8786337147_wp,  0.8184929719_wp, & !53-56
      &  0.8272905671_wp,  0.1792896594_wp,  0.2333308237_wp,  0.3160232084_wp, & !57-60
      &  0.2048800101_wp,  0.6103831629_wp,  0.2847769346_wp,  0.3365440133_wp, & !61-64
      &  0.4473717293_wp,  0.3302940420_wp,  0.3465382047_wp,  0.3117739683_wp, & !65-68
      &  0.3320431148_wp,  0.2821686524_wp,  0.2700901299_wp,  0.7547418385_wp, & !69-72
      &  0.6671709344_wp,  1.2829896935_wp,  1.5071299849_wp,  1.0844648293_wp, & !73-76
      &  1.3253676523_wp,  2.1863745692_wp,  1.7200032152_wp,  1.0036728849_wp, & !77-80
      &  0.7714939233_wp,  1.0695806778_wp,  1.3167890037_wp,  0.9040047383_wp, & !81-84
      &  1.2378987054_wp,  0.3214306997_wp,  0.2322560689_wp,  0.6547977477_wp, & !85-88
      &  0.5164282791_wp,  0.4025001070_wp,  0.3624847114_wp,  0.4383068517_wp, & !89-92
      &  0.4782382905_wp,  0.4877346759_wp,  0.5713826416_wp,  0.5899939524_wp, & !93-96
      &  0.6286367467_wp,  0.5236962126_wp,  0.6522982637_wp,  0.6940382810_wp, & !97-100
      &  0.4244551058_wp,  0.3569789853_wp,  0.5908324555_wp] !101-103

   !> Element-specific covalent radii for the CN for the EEQ_BC charges.
   real(wp), parameter :: eeqbc_cov_radii(max_elem) = 0.5_wp*[&
      &  0.6429899185_wp,  0.7291612708_wp,  2.5202214796_wp,  2.2184411427_wp, & !1-4
      &  2.1072579060_wp,  2.1978094892_wp,  2.4070801098_wp,  2.4073999877_wp, & !5-8
      &  2.2875231881_wp,  1.6311402061_wp,  3.1771542577_wp,  3.0074402025_wp, & !9-12
      &  2.3697795763_wp,  2.9498019121_wp,  3.2281425082_wp,  3.3769351202_wp, & !13-16
      &  3.1425408157_wp,  1.5111948734_wp,  3.1923508047_wp,  3.1729127730_wp, & !17-20
      &  2.8894792144_wp,  2.9799903432_wp,  2.7166595480_wp,  2.7032097879_wp, & !21-24
      &  2.4640716466_wp,  3.0487196769_wp,  3.1415235617_wp,  3.0632129440_wp, & !25-28
      &  2.3437338422_wp,  2.8223903756_wp,  3.0673440934_wp,  3.2673003020_wp, & !29-32
      &  3.4446562432_wp,  3.5308120404_wp,  3.5651158793_wp,  2.3735083821_wp, & !33-36
      &  3.5739640685_wp,  3.3777302404_wp,  2.9943323822_wp,  3.3564146158_wp, & !37-40
      &  3.1851646496_wp,  3.2113979549_wp,  3.0056059341_wp,  3.3481004957_wp, & !41-44
      &  3.3048651615_wp,  3.2740005664_wp,  3.0837075099_wp,  3.2584736917_wp, & !45-48
      &  3.3214672536_wp,  3.7304951804_wp,  3.6004609601_wp,  3.8608037236_wp, & !49-52
      &  3.9141725828_wp,  2.6862414414_wp,  3.8077629678_wp,  3.8010046557_wp, & !53-56
      &  2.7806681846_wp,  2.6467597727_wp,  3.6526201038_wp,  3.9132683142_wp, & !57-60
      &  3.8239010987_wp,  2.9373789663_wp,  3.4787640374_wp,  3.3993270571_wp, & !61-64
      &  3.6017484887_wp,  3.6759255923_wp,  3.6460860950_wp,  3.8275261532_wp, & !65-68
      &  3.7282819662_wp,  3.6475356348_wp,  3.4415844962_wp,  3.1763798872_wp, & !69-72
      &  2.9069418994_wp,  3.1465466911_wp,  2.8525944634_wp,  3.3578656946_wp, & !73-76
      &  3.4024517035_wp,  3.5417424786_wp,  3.1202031706_wp,  3.6892972439_wp, & !77-80
      &  3.5049543259_wp,  3.9944371767_wp,  3.4954910114_wp,  3.8425468344_wp, & !81-84
      &  4.1301856276_wp,  2.9340156007_wp,  3.9620682706_wp,  3.5310183521_wp, & !85-88
      &  1.7855486856_wp,  1.1206587951_wp,  1.6449496462_wp,  2.7089924333_wp, & !89-92
      &  2.3890886476_wp,  3.0906012363_wp,  2.8396513736_wp,  2.8300586548_wp, & !93-96
      &  2.4717129593_wp,  2.7263021764_wp,  3.0489604535_wp,  1.3684588197_wp, & !97-100
      &  2.6129438461_wp,  2.7433186093_wp,  2.6699912181_wp] !101-103

   !> Element-specific averaged coordination number over the fitset for the EEQ_BC charges.
   real(wp), parameter :: eeqbc_avg_cn(max_elem) = [&
      &  0.3921100000_wp, 0.0810600000_wp, 0.9910100000_wp, 0.7499500000_wp, & !1-4
      &  1.1543700000_wp, 1.6691400000_wp, 1.4250300000_wp, 0.8718100000_wp, & !5-8
      &  0.6334000000_wp, 0.0876700000_wp, 0.8740600000_wp, 0.8754800000_wp, & !9-12
      &  1.2147200000_wp, 1.1335000000_wp, 1.6890600000_wp, 1.0221600000_wp, & !13-16
      &  0.5386400000_wp, 0.0827800000_wp, 1.4096300000_wp, 1.1954700000_wp, & !17-20
      &  1.5142100000_wp, 1.7892000000_wp, 2.0646100000_wp, 1.6905600000_wp, & !21-24
      &  1.6563700000_wp, 1.5128400000_wp, 1.3179000000_wp, 0.9749800000_wp, & !25-28
      &  0.5334600000_wp, 0.6585000000_wp, 0.9696500000_wp, 1.0083100000_wp, & !29-32
      &  1.0871000000_wp, 0.8222200000_wp, 0.5449300000_wp, 0.1647100000_wp, & !33-36
      &  1.2490800000_wp, 1.2198700000_wp, 1.5657400000_wp, 1.8697600000_wp, & !37-40
      &  1.8947900000_wp, 1.7085000000_wp, 1.5521300000_wp, 1.4903300000_wp, & !41-44
      &  1.3177400000_wp, 0.6991700000_wp, 0.5528200000_wp, 0.6642200000_wp, & !45-48
      &  0.9069800000_wp, 1.0976200000_wp, 1.2183000000_wp, 0.7321900000_wp, & !49-52
      &  0.5498700000_wp, 0.2467100000_wp, 1.5680600000_wp, 1.1677300000_wp, & !53-56
      &  1.6642500000_wp, 1.6032600000_wp, 1.6032600000_wp, 1.6032600000_wp, & !57-60
      &  1.6032600000_wp, 1.6032600000_wp, 1.6032600000_wp, 1.6032600000_wp, & !61-64
      &  1.6032600000_wp, 1.6032600000_wp, 1.6032600000_wp, 1.6032600000_wp, & !65-68
      &  1.6032600000_wp, 1.6032600000_wp, 1.6032600000_wp, 1.8191000000_wp, & !69-72
      &  1.8175100000_wp, 1.6802300000_wp, 1.5224100000_wp, 1.4602600000_wp, & !73-76
      &  1.1110400000_wp, 0.9102600000_wp, 0.5218000000_wp, 1.4895900000_wp, & !77-80
      &  0.8441800000_wp, 0.9426900000_wp, 1.5171900000_wp, 0.7287100000_wp, & !81-84
      &  0.5137000000_wp, 0.2678200000_wp, 1.2122500000_wp, 1.5797100000_wp, & !85-88
      &  1.7549800000_wp, 1.7549800000_wp, 1.7549800000_wp, 1.7549800000_wp, & !89-92
      &  1.7549800000_wp, 1.7549800000_wp, 1.7549800000_wp, 1.7549800000_wp, & !93-96
      &  1.7549800000_wp, 1.7549800000_wp, 1.7549800000_wp, 1.7549800000_wp, & !97-100
      &  1.7549800000_wp, 1.7549800000_wp, 1.7549800000_wp] !101-103

   !> Element-specific scaling of the pairwise van-der-Waals radii.
   real(wp), parameter :: eeqbc_rvdw_scale(max_elem) = [&
      &  1.0000000000_wp,  1.0000000000_wp,  1.0500000000_wp,  1.0000000000_wp, & !1-4
      &  1.0500000000_wp,  1.0000000000_wp,  1.0500000000_wp,  1.0600000000_wp, & !5-8
      &  1.1700000000_wp,  0.7500000000_wp,  1.0800000000_wp,  1.1300000000_wp, & !9-12
      &  1.1300000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !13-16
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !17-20
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !21-24
      &  1.0000000000_wp,  1.1000000000_wp,  1.1000000000_wp,  1.1000000000_wp, & !25-28
      &  1.0000000000_wp,  1.1500000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !29-32
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !33-36
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !37-40
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !41-44
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.1500000000_wp, & !45-48
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !49-52
      &  1.0000000000_wp,  1.0000000000_wp,  0.5500000000_wp,  1.0000000000_wp, & !53-56
      &  1.2500000000_wp,  1.2500000000_wp,  1.2000000000_wp,  1.2000000000_wp, & !57-60
      &  1.2000000000_wp,  1.0000000000_wp,  1.2000000000_wp,  1.2000000000_wp, & !61-64
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !65-68
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !69-72
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !73-76
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.1500000000_wp, & !77-80
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !81-84
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !85-88
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !89-92
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !93-96
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !97-100
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp] !101-103

contains

!> Get electronegativity for species with a given symbol
elemental function get_eeqbc_chi_sym(symbol) result(chi)

   !> Element symbol
   character(len=*), intent(in) :: symbol

   !> electronegativity
   real(wp) :: chi

   chi = get_eeqbc_chi(to_number(symbol))

end function get_eeqbc_chi_sym

!> Get electronegativity for species with a given atomic number
elemental function get_eeqbc_chi_num(number) result(chi)

   !> Atomic number
   integer, intent(in) :: number

   !> electronegativity
   real(wp) :: chi

   if (number > 0 .and. number <= size(eeqbc_chi, dim=1)) then
      chi = eeqbc_chi(number)
   else
      chi = -1.0_wp
   end if

end function get_eeqbc_chi_num

!> Get hardness for species with a given symbol
elemental function get_eeqbc_eta_sym(symbol) result(eta)

   !> Element symbol
   character(len=*), intent(in) :: symbol

   !> hardness
   real(wp) :: eta

   eta = get_eeqbc_eta(to_number(symbol))

end function get_eeqbc_eta_sym

!> Get hardness for species with a given atomic number
elemental function get_eeqbc_eta_num(number) result(eta)

   !> Atomic number
   integer, intent(in) :: number

   !> hardness
   real(wp) :: eta

   if (number > 0 .and. number <= size(eeqbc_eta, dim=1)) then
      eta = eeqbc_eta(number)
   else
      eta = -1.0_wp
   end if

end function get_eeqbc_eta_num

!> Get charge width for species with a given symbol
elemental function get_eeqbc_rad_sym(symbol) result(rad)

   !> Element symbol
   character(len=*), intent(in) :: symbol

   !> charge width
   real(wp) :: rad

   rad = get_eeqbc_rad(to_number(symbol))

end function get_eeqbc_rad_sym

!> Get charge width for species with a given atomic number
elemental function get_eeqbc_rad_num(number) result(rad)

   !> Atomic number
   integer, intent(in) :: number

   !> charge width
   real(wp) :: rad

   if (number > 0 .and. number <= size(eeqbc_rad, dim=1)) then
      rad = eeqbc_rad(number)
   else
      rad = -1.0_wp
   end if

end function get_eeqbc_rad_num

!> Get CN scaling of the electronegativity for species with a given symbol
elemental function get_eeqbc_kcnchi_sym(symbol) result(kcnchi)

   !> Element symbol
   character(len=*), intent(in) :: symbol

   !> CN scaling of EN
   real(wp) :: kcnchi

   kcnchi = get_eeqbc_kcnchi(to_number(symbol))

end function get_eeqbc_kcnchi_sym

!> Get CN scaling of the electronegativity for species with a given atomic number
elemental function get_eeqbc_kcnchi_num(number) result(kcnchi)

   !> Atomic number
   integer, intent(in) :: number

   !> CN scaling of EN
   real(wp) :: kcnchi

   if (number > 0 .and. number <= size(eeqbc_kcnchi, dim=1)) then
      kcnchi = eeqbc_kcnchi(number)
   else
      kcnchi = -1.0_wp
   end if

end function get_eeqbc_kcnchi_num

!> Get local q scaling of the electronegativity for species with a given symbol
elemental function get_eeqbc_kqchi_sym(symbol) result(kqchi)

   !> Element symbol
   character(len=*), intent(in) :: symbol

   !> local q scaling of EN
   real(wp) :: kqchi

   kqchi = get_eeqbc_kqchi(to_number(symbol))

end function get_eeqbc_kqchi_sym

!> Get local q scaling of the electronegativity for species with a given atomic number
elemental function get_eeqbc_kqchi_num(number) result(kqchi)

   !> Atomic number
   integer, intent(in) :: number

   !> local q scaling of EN
   real(wp) :: kqchi

   if (number > 0 .and. number <= size(eeqbc_kqchi, dim=1)) then
      kqchi = eeqbc_kqchi(number)
   else
      kqchi = -1.0_wp
   end if

end function get_eeqbc_kqchi_num

!> Get local q scaling of the chemical hardness for species with a given symbol
elemental function get_eeqbc_kqeta_sym(symbol) result(kqeta)

   !> Element symbol
   character(len=*), intent(in) :: symbol

   !> local q scaling of hardness
   real(wp) :: kqeta

   kqeta = get_eeqbc_kqeta(to_number(symbol))

end function get_eeqbc_kqeta_sym

!> Get local q scaling of the chemical hardness for species with a given atomic number
elemental function get_eeqbc_kqeta_num(number) result(kqeta)

   !> Atomic number
   integer, intent(in) :: number

   !> local q scaling of hardness
   real(wp) :: kqeta

   if (number > 0 .and. number <= size(eeqbc_kqeta, dim=1)) then
      kqeta = eeqbc_kqeta(number)
   else
      kqeta = -1.0_wp
   end if

end function get_eeqbc_kqeta_num

!> Get bond capacitance for species with a given symbol
elemental function get_eeqbc_cap_sym(symbol) result(cap)

   !> Element symbol
   character(len=*), intent(in) :: symbol

   !> bond capacitance
   real(wp) :: cap

   cap = get_eeqbc_cap(to_number(symbol))

end function get_eeqbc_cap_sym

!> Get bond capacitance for species with a given atomic number
elemental function get_eeqbc_cap_num(number) result(cap)

   !> Atomic number
   integer, intent(in) :: number

   !> bond capacitance
   real(wp) :: cap

   if (number > 0 .and. number <= size(eeqbc_cap, dim=1)) then
      cap = eeqbc_cap(number)
   else
      cap = -1.0_wp
   end if

end function get_eeqbc_cap_num

!> Get covalent radius for species with a given symbol
elemental function get_eeqbc_cov_radii_sym(symbol) result(rcov)

   !> Element symbol
   character(len=*), intent(in) :: symbol

   !> covalent radius
   real(wp) :: rcov

   rcov = get_eeqbc_cov_radii(to_number(symbol))

end function get_eeqbc_cov_radii_sym

!> Get covalent radius for species with a given atomic number
elemental function get_eeqbc_cov_radii_num(number) result(rcov)

   !> Atomic number
   integer, intent(in) :: number

   !> covalent radius
   real(wp) :: rcov

   if (number > 0 .and. number <= size(eeqbc_cov_radii, dim=1)) then
      rcov = eeqbc_cov_radii(number)
   else
      rcov = -1.0_wp
   end if

end function get_eeqbc_cov_radii_num

!> Get average CN for species with a given symbol
elemental function get_eeqbc_avg_cn_sym(symbol) result(avg_cn)

   !> Element symbol
   character(len=*), intent(in) :: symbol

   !> average CN
   real(wp) :: avg_cn

   avg_cn = get_eeqbc_avg_cn(to_number(symbol))

end function get_eeqbc_avg_cn_sym

!> Get average CN for species with a given atomic number
elemental function get_eeqbc_avg_cn_num(number) result(avg_cn)

   !> Atomic number
   integer, intent(in) :: number

   !> average CN
   real(wp) :: avg_cn

   if (number > 0 .and. number <= size(eeqbc_avg_cn, dim=1)) then
      avg_cn = eeqbc_avg_cn(number)
   else
      avg_cn = -1.0_wp
   end if

end function get_eeqbc_avg_cn_num

!> Get scaling for pairwise van-der-Waals radius for species with a given symbol
elemental function get_eeqbc_rvdw_scale_sym(symbol) result(rvdw_scale)

   !> Element symbol
   character(len=*), intent(in) :: symbol

   !> scaling for van-der-Waals radius
   real(wp) :: rvdw_scale

   rvdw_scale = get_eeqbc_rvdw_scale(to_number(symbol))

end function get_eeqbc_rvdw_scale_sym

!> Get scaling for pairwise van-der-Waals radius for species with a given atomic number
elemental function get_eeqbc_rvdw_scale_num(number) result(rvdw_scale)

   !> Atomic number
   integer, intent(in) :: number

   !> scaling for van-der-Waals radius
   real(wp) :: rvdw_scale

   if (number > 0 .and. number <= size(eeqbc_rvdw_scale, dim=1)) then
      rvdw_scale = eeqbc_rvdw_scale(number)
   else
      rvdw_scale = -1.0_wp
   end if

end function get_eeqbc_rvdw_scale_num

end module multicharge_param_eeqbc2025
