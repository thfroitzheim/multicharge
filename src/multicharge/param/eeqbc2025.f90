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
!> The parametrization was revised due to minor inconsistencies mainly
!> for the noble gases. In this course the charge dependency of the 
!> hardness value could be eliminated to reduce the element-wise parameters
!> to seven. The new parametrization should generally be preferred.
!> The original parametrization can be used in multicharge v0.5.0. 
!> Details on the modifications are published in
!>
!> Benedikt Bädorf, ..., Thomas Froitzheim, Stefan Grimme, ...
module multicharge_param_eeqbc2025
   use mctc_env, only: wp
   use mctc_io_symbols, only: to_number
   implicit none
   private

   public :: get_eeqbc_chi, get_eeqbc_eta, get_eeqbc_rad, get_eeqbc_kcnchi, &
      & get_eeqbc_kqchi, get_eeqbc_cap, get_eeqbc_cov_radii, get_eeqbc_avg_cn, &
      & get_eeqbc_rvdw_scale

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
      &  1.5191499066_wp, -0.9209695912_wp, -1.0984929066_wp,  0.0608063848_wp, & !1-4
      &  0.5774005116_wp,  1.4041145025_wp,  1.7983973228_wp,  1.9680353460_wp, & !5-8
      &  1.6512764169_wp, -1.0631757735_wp, -1.3422793213_wp, -0.6069600018_wp, & !9-12
      & -0.1112475104_wp,  0.2856667465_wp,  1.2000252650_wp,  1.6046422023_wp, & !13-16
      &  1.5326006552_wp, -0.8963941419_wp, -1.8383037641_wp, -0.8259272000_wp, & !17-20
      & -0.7497398690_wp, -0.4091948277_wp, -0.2864952025_wp, -0.7283757774_wp, & !21-24
      & -0.1618705534_wp,  0.0315198954_wp,  0.3326673500_wp,  0.3090962564_wp, & !25-28
      & -0.1926813854_wp, -0.3312483278_wp, -0.2588647514_wp,  0.2509840410_wp, & !29-32
      &  0.6526405847_wp,  1.5505499977_wp,  1.2044205334_wp, -1.1249298043_wp, & !33-36
      & -2.0407952931_wp, -0.7988994644_wp, -0.6350541351_wp, -0.5258877375_wp, & !37-40
      & -0.4357000229_wp, -0.3619492361_wp, -0.0736777235_wp, -0.1925290646_wp, & !41-44
      &  0.1058976020_wp, -0.1681250797_wp, -0.2272611425_wp, -0.2560662966_wp, & !45-48
      & -0.4152287683_wp,  0.1561360304_wp,  0.3136925937_wp,  1.3256154015_wp, & !49-52
      &  0.9170416132_wp, -0.8922326611_wp, -2.2291629633_wp, -1.0518337127_wp, & !53-56
      & -0.9906000446_wp, -1.2213400352_wp, -2.0454701071_wp, -1.9476145462_wp, & !57-60
      & -2.0620178988_wp, -1.7367188677_wp, -1.6452965843_wp, -1.3172436591_wp, & !61-64
      & -1.6504494236_wp, -1.8910989723_wp, -1.8660574339_wp, -1.8017150448_wp, & !65-68
      & -1.8238057774_wp, -0.8732520767_wp, -0.6528872798_wp, -0.5404585056_wp, & !69-72
      & -0.3079337519_wp,  0.0196214635_wp,  0.1247887764_wp,  0.2056717454_wp, & !73-76
      &  0.3021607718_wp,  0.3912207224_wp,  0.2789923468_wp, -0.3331306716_wp, & !77-80
      & -0.7700269731_wp, -0.2660374186_wp,  0.2392071619_wp,  0.8817638893_wp, & !81-84
      &  0.9358154209_wp, -0.6042651825_wp, -1.9569712252_wp, -0.9546435603_wp, & !85-88
      & -0.6774106112_wp, -0.5004270055_wp, -0.4677650154_wp, -0.7509765459_wp, & !89-92
      & -0.2654613531_wp, -0.7397485542_wp, -0.5125040933_wp, -0.5002606895_wp, & !93-96
      & -0.4357129799_wp, -0.7942980087_wp, -0.6116598690_wp, -0.5795357528_wp, & !97-100
      & -0.5927970666_wp, -0.5224677395_wp, -0.7921783756_wp] !101-103

   !> Element-specific chemical hardnesses for the EEQ_BC charges.
   real(wp), parameter :: eeqbc_eta(max_elem) = [&
      &  9.0582061511_wp,  0.8041649428_wp,  0.7583095441_wp,  1.8078533790_wp, & !1-4
      & -1.7187777926_wp, -3.2271175403_wp, -1.7746769474_wp, -0.5638967716_wp, & !5-8
      & -2.5977658763_wp, -0.0787334753_wp, -1.2402083366_wp, -0.8258891523_wp, & !9-12
      & -1.2876451559_wp, -2.7051676357_wp, -8.5282941664_wp, -2.0018869403_wp, & !13-16
      &  3.7262890011_wp,  1.4760185605_wp,  1.1718388905_wp, -0.0160361511_wp, & !17-20
      & -0.8017191239_wp, -0.6507438424_wp, -0.4314604838_wp,  0.8718524059_wp, & !21-24
      &  1.3348986867_wp, -1.2361272649_wp, -1.3286739672_wp, -0.6029795268_wp, & !25-28
      &  0.8485394404_wp,  0.2016004558_wp,  1.7374086049_wp, -1.0369979205_wp, & !29-32
      & -1.1353279286_wp,  1.2256462285_wp,  2.2551606799_wp,  0.9921918377_wp, & !33-36
      &  1.8027877838_wp, -1.5593199936_wp, -0.6678216168_wp, -0.4885755038_wp, & !37-40
      & -1.2908022726_wp,  0.2370153718_wp,  1.3177903990_wp,  1.0749577623_wp, & !41-44
      &  0.0501671847_wp,  0.6957545544_wp,  2.2063519566_wp, -0.1835965684_wp, & !45-48
      &  1.0107253626_wp, -1.1294660856_wp, -1.4629726062_wp,  0.5799218982_wp, & !49-52
      &  1.5313833959_wp,  2.7665693682_wp,  1.3941409467_wp, -0.8230792968_wp, & !53-56
      & -0.8462414258_wp, -0.8956236640_wp, -0.4520112608_wp, -0.6870875169_wp, & !57-60
      & -0.8148495356_wp, -0.9915398260_wp, -1.1144089642_wp, -1.3415909696_wp, & !61-64
      & -1.2743913536_wp, -0.8475833340_wp, -0.8235772327_wp, -0.6134570702_wp, & !65-68
      & -0.2057688045_wp, -0.1970843227_wp,  0.2283553963_wp, -1.1368202968_wp, & !69-72
      & -0.6488618911_wp,  0.5035417313_wp,  1.2125579931_wp,  0.9126551240_wp, & !73-76
      &  1.0837059704_wp,  1.1920324217_wp,  0.5069340207_wp,  2.3407278834_wp, & !77-80
      &  1.7075772993_wp, -0.2080852805_wp, -5.2344801158_wp, -0.8192834449_wp, & !81-84
      &  1.3743911092_wp,  4.0909791895_wp,  0.5932976353_wp, -1.1924443432_wp, & !85-88
      & -1.0189590540_wp, -0.4109049288_wp, -0.3578600129_wp,  0.0575588392_wp, & !89-92
      & -0.0996545934_wp, -0.2423792004_wp, -0.2815998574_wp, -0.2893275429_wp, & !93-96
      & -0.1463647915_wp, -0.1939367318_wp, -0.3900228407_wp, -0.7418098647_wp, & !97-100
      & -0.8421135587_wp, -0.9730365716_wp, -1.7585136827_wp] !101-103

   !> Element-specific charge widths for the EEQ_BC charges.
   real(wp), parameter :: eeqbc_rad(max_elem) = [&
      &  0.9839410084_wp,  0.0217134247_wp,  0.1638324678_wp,  0.3609269021_wp, & !1-4
      &  0.1579359051_wp,  0.1365739315_wp,  0.2208151103_wp,  0.2531042155_wp, & !5-8
      &  0.1211682446_wp,  0.0338374633_wp,  0.1414484758_wp,  0.2099057654_wp, & !9-12
      &  0.1755175287_wp,  0.1632670964_wp,  0.0750710068_wp,  0.1924626897_wp, & !13-16
      &  0.4209371540_wp,  0.1093954534_wp,  0.2274606824_wp,  0.3401629543_wp, & !17-20
      &  0.2597845352_wp,  0.2754819656_wp,  0.2763612141_wp,  0.4061457983_wp, & !21-24
      &  0.6687264806_wp,  0.2690020498_wp,  0.2487876389_wp,  0.2547007418_wp, & !25-28
      &  0.3262174130_wp,  0.2849786989_wp,  0.5019989767_wp,  0.2569084936_wp, & !29-32
      &  0.2865182174_wp,  0.4788098053_wp,  0.5422437853_wp,  0.1472678593_wp, & !33-36
      &  0.2814149796_wp,  0.2167931645_wp,  0.3137623033_wp,  0.3139276651_wp, & !37-40
      &  0.2444021175_wp,  0.4414022629_wp,  0.8627055315_wp,  0.4767406420_wp, & !41-44
      &  0.3833942423_wp,  0.3197155801_wp,  0.6253483607_wp,  0.2330676679_wp, & !45-48
      &  0.3584219499_wp,  0.2624412731_wp,  0.2224717030_wp,  0.4490757329_wp, & !49-52
      &  0.5941559606_wp,  0.3687821430_wp,  0.2791474685_wp,  0.2762801484_wp, & !53-56
      &  0.2913437480_wp,  0.2774555696_wp,  0.2289366469_wp,  0.2292336683_wp, & !57-60
      &  0.2150517934_wp,  0.2088650053_wp,  0.2017240195_wp,  0.2016399205_wp, & !61-64
      &  0.2025157285_wp,  0.1964457976_wp,  0.2167650199_wp,  0.2420319819_wp, & !65-68
      &  0.2651957474_wp,  0.3466690461_wp,  0.4281820135_wp,  0.2806381217_wp, & !69-72
      &  0.3046837831_wp,  0.5462094680_wp,  0.6525123701_wp,  0.6220281066_wp, & !73-76
      &  0.6383362450_wp,  0.5228488673_wp,  0.5065444496_wp,  0.5190734756_wp, & !77-80
      &  0.4316437065_wp,  0.3644215444_wp,  0.1205035672_wp,  0.3501864560_wp, & !81-84
      &  0.6458601028_wp,  1.0952923287_wp,  0.2330814226_wp,  0.2478659987_wp, & !85-88
      &  0.2950302783_wp,  0.3353614900_wp,  0.3914166269_wp,  0.4151244844_wp, & !89-92
      &  0.4514498186_wp,  0.3581013419_wp,  0.3704320608_wp,  0.3643862940_wp, & !93-96
      &  0.4340008100_wp,  0.3707125334_wp,  0.3616211901_wp,  0.3237746386_wp, & !97-100
      &  0.2891314661_wp,  0.2419831064_wp,  0.2074893721_wp] !101-103

   !> Element-specific CN scaling of the electronegativity for the EEQ_BC charges.
   real(wp), parameter :: eeqbc_kcnchi(max_elem) = [&
      &  4.4968302538_wp,  6.7617987418_wp, -1.3171759988_wp, -0.4790274326_wp, & !1-4
      & -0.2115016269_wp,  0.8208173413_wp,  1.1901660433_wp,  1.6642055948_wp, & !5-8
      &  2.8570548551_wp, 13.1259343303_wp, -1.0261716631_wp, -0.5514568804_wp, & !9-12
      &  0.3010242004_wp,  0.2283503050_wp,  1.4665587614_wp,  1.8958844849_wp, & !13-16
      &  3.9568521008_wp, 17.0808375779_wp, -2.2896763736_wp, -2.1846866044_wp, & !17-20
      & -1.5796385849_wp, -1.4443577316_wp, -1.2453088832_wp, -0.7963902696_wp, & !21-24
      & -1.8272432717_wp, -1.2595036774_wp, -1.0601372114_wp, -1.3797547956_wp, & !25-28
      &  0.1858351472_wp,  0.0601068231_wp, -0.1293734233_wp,  0.1283979831_wp, & !29-32
      &  0.8055642048_wp,  2.0705415705_wp,  3.0312084391_wp,  7.9819897255_wp, & !33-36
      & -2.4305200236_wp, -2.0967233879_wp, -1.8823894915_wp, -1.4124569058_wp, & !37-40
      & -0.4384883546_wp, -0.4101768976_wp, -1.0469185277_wp, -0.5218530820_wp, & !41-44
      & -0.2986745225_wp,  1.1273402854_wp,  0.0762684616_wp,  0.7056137470_wp, & !45-48
      & -0.1607980287_wp,  0.2279659476_wp,  0.5531340654_wp,  1.7626542489_wp, & !49-52
      &  2.4347180522_wp,  3.8452238473_wp, -3.4375912689_wp, -2.4651282340_wp, & !53-56
      & -2.1518205531_wp, -2.3127120902_wp, -2.5610116467_wp, -2.6748326091_wp, & !57-60
      & -2.8152071856_wp, -2.1343549190_wp, -1.4577812067_wp, -2.3523465836_wp, & !61-64
      & -1.9352942923_wp, -1.2660821881_wp, -2.0963638687_wp, -2.0933816272_wp, & !65-68
      & -1.9439381912_wp, -1.0249686067_wp, -1.8397722797_wp, -1.6617519771_wp, & !69-72
      & -0.9087050758_wp, -0.3116399072_wp, -0.8195612347_wp, -0.6722865351_wp, & !73-76
      & -0.6062531950_wp,  0.3092767176_wp,  0.9677260442_wp,  0.4771060067_wp, & !77-80
      & -0.6800762246_wp,  0.0211766650_wp,  0.4178613504_wp,  1.2759950786_wp, & !81-84
      &  2.2409755176_wp,  3.0411789217_wp, -3.4803010946_wp, -2.3237228487_wp, & !85-88
      & -2.6141694563_wp, -3.6612802249_wp, -3.1173088793_wp, -2.7965359657_wp, & !89-92
      & -3.3254939886_wp, -2.6532628418_wp, -2.9391713100_wp, -3.1945879627_wp, & !93-96
      & -2.4167872826_wp, -2.5102597158_wp, -1.8356056600_wp, -2.4744474703_wp, & !97-100
      & -1.9648309065_wp, -1.6882157971_wp, -2.8688650895_wp] !101-103

   !> Element-specific local q scaling of the electronegativity for the EEQ_BC charges.
   real(wp), parameter :: eeqbc_kqchi(max_elem) = [&
      &  1.3303531125_wp, -1.7242134420_wp,  5.7705977770_wp,  4.6626931398_wp, & !1-4
      &  4.4534579986_wp,  2.9350024757_wp,  1.6081200500_wp,  1.7424233175_wp, & !5-8
      &  1.5032457384_wp, -0.8765009932_wp,  5.2107736875_wp,  4.4831159022_wp, & !9-12
      &  4.4536135519_wp,  3.7680938785_wp,  3.6217987809_wp,  2.6922230774_wp, & !13-16
      &  3.4649839240_wp,  0.2468249576_wp,  4.9358350765_wp,  4.8526712772_wp, & !17-20
      &  4.2658096082_wp,  4.5502257175_wp,  4.2659758527_wp,  3.2196366336_wp, & !21-24
      &  3.9959795437_wp,  3.1846380268_wp,  3.7488609058_wp,  4.5265085862_wp, & !25-28
      &  3.8508249527_wp,  3.6743418339_wp,  4.2227509728_wp,  3.9622581355_wp, & !29-32
      &  2.8188203533_wp,  2.8442925258_wp,  2.9344345980_wp,  0.0199520155_wp, & !33-36
      &  5.8780737052_wp,  5.3357767737_wp,  4.6056885121_wp,  4.3653701096_wp, & !37-40
      &  3.8271654393_wp,  3.9355522733_wp,  3.5514542225_wp,  4.2538691630_wp, & !41-44
      &  4.2896075117_wp,  3.1982871310_wp,  3.5906891069_wp,  3.8838109900_wp, & !45-48
      &  4.5835267808_wp,  3.8033548642_wp,  3.4869440792_wp,  2.9487263998_wp, & !49-52
      &  2.9465360157_wp,  1.1512965161_wp,  4.8916809590_wp,  5.3007095014_wp, & !53-56
      &  4.6315507022_wp,  4.3062639119_wp,  5.7896863535_wp,  5.4854825010_wp, & !57-60
      &  5.6893192836_wp,  5.2603769526_wp,  4.7881893218_wp,  4.9982455488_wp, & !61-64
      &  4.9072809771_wp,  4.8883176185_wp,  5.1750378792_wp,  4.5929043252_wp, & !65-68
      &  4.7474842549_wp,  3.8133530664_wp,  4.2977685770_wp,  4.3139800288_wp, & !69-72
      &  3.9504487197_wp,  4.0309384009_wp,  3.6221491705_wp,  4.1656314720_wp, & !73-76
      &  3.8192046831_wp,  3.9236746345_wp,  3.4017377294_wp,  3.3528043961_wp, & !77-80
      &  5.0413392645_wp,  4.2940947806_wp,  3.8134931081_wp,  3.0080683902_wp, & !81-84
      &  3.1403625651_wp,  1.4527821515_wp,  5.8643808717_wp,  4.8332445767_wp, & !85-88
      &  4.9109443197_wp,  5.9854518669_wp,  6.1519484325_wp,  5.2292772742_wp, & !89-92
      &  6.3678703643_wp,  5.4082399086_wp,  5.8591145923_wp,  5.6866076813_wp, & !93-96
      &  5.7601172907_wp,  5.6654685712_wp,  5.4101462804_wp,  5.2220228661_wp, & !97-100
      &  4.8570597562_wp,  5.4523814609_wp,  5.5738399131_wp] !101-103

   !> Element-specific bond capacitance for the EEQ_BC charges.
   real(wp), parameter :: eeqbc_cap(max_elem) = [&
      &  2.5216276196_wp,  0.5752012095_wp,  1.2369573243_wp,  1.4568108519_wp, & !1-4
      &  2.2781791024_wp,  6.1415957539_wp,  5.9564663025_wp,  7.3478606982_wp, & !5-8
      &  6.6896494671_wp,  2.3995369155_wp,  2.0494434593_wp,  1.3634808882_wp, & !9-12
      &  0.6764709418_wp,  1.0553930553_wp,  0.5838983269_wp,  1.0698694628_wp, & !13-16
      &  2.7228310916_wp,  0.8880788013_wp,  2.1036490840_wp,  2.1745284033_wp, & !17-20
      &  2.1220263472_wp,  3.9600967119_wp,  2.1951678870_wp,  1.8083448239_wp, & !21-24
      &  3.7542002877_wp,  2.6587878078_wp,  4.9166994447_wp,  5.8478097753_wp, & !25-28
      &  3.5216040648_wp,  2.0771852618_wp,  1.1816592215_wp,  1.9425542033_wp, & !29-32
      &  1.1886785638_wp,  0.9907398847_wp,  3.2313908893_wp,  1.7834502553_wp, & !33-36
      &  1.9923237905_wp,  3.5305347643_wp,  3.9807682911_wp,  2.3004053992_wp, & !37-40
      &  1.8006878044_wp,  2.6174128574_wp,  2.8385492577_wp,  5.9511023128_wp, & !41-44
      &  7.3219068265_wp,  2.5946677227_wp,  2.4657310186_wp,  1.1521860214_wp, & !45-48
      &  1.6359304350_wp,  1.9200311523_wp,  0.9497476961_wp,  1.1954053101_wp, & !49-52
      &  1.9587922570_wp,  3.0380808310_wp,  2.1733977084_wp,  1.5382385588_wp, & !53-56
      &  1.8488981879_wp,  1.2317159272_wp,  0.5400245361_wp,  0.6863530982_wp, & !57-60
      &  0.6226898464_wp,  0.6243838696_wp,  0.5391835837_wp,  0.9920311462_wp, & !61-64
      &  1.2182253771_wp,  0.6487679115_wp,  0.8954406813_wp,  1.1176587704_wp, & !65-68
      &  0.8749487933_wp,  2.4706672088_wp,  4.0623677775_wp,  3.3940889516_wp, & !69-72
      &  0.8810673853_wp,  3.7121712038_wp,  3.3579787640_wp,  2.5318266269_wp, & !73-76
      &  1.5700999815_wp,  3.4359384954_wp,  9.4616163161_wp,  3.0663498188_wp, & !77-80
      &  0.8326138553_wp,  2.1710985520_wp,  3.2739620216_wp,  1.0430012664_wp, & !81-84
      &  6.0946540971_wp,  2.3157509428_wp,  2.3427591940_wp,  2.6297147306_wp, & !85-88
      &  1.3429457496_wp,  0.5516647402_wp,  3.3187434315_wp,  2.4653267014_wp, & !89-92
      &  5.4148346306_wp,  3.9178141298_wp,  6.0473116069_wp,  8.7626930932_wp, & !93-96
      &  6.1503056318_wp,  2.9568162565_wp,  4.5100302745_wp,  2.8271998289_wp, & !97-100
      &  1.4311502993_wp,  1.9220734974_wp,  1.6493745255_wp] !101-103

   !> Element-specific covalent radii for the CN for the EEQ_BC charges.
   real(wp), parameter :: eeqbc_cov_radii(max_elem) = 0.5_wp*[&
      &  1.0873678902_wp,  0.0045628280_wp,  2.8385414023_wp,  2.2369359793_wp, & !1-4
      &  2.2631432568_wp,  2.5556464299_wp,  2.6528219471_wp,  2.5471166478_wp, & !5-8
      &  2.0970520036_wp,  1.1527679853_wp,  3.9222564151_wp,  3.6628112720_wp, & !9-12
      &  3.1200757440_wp,  3.2311571633_wp,  3.3714412240_wp,  3.4966508157_wp, & !13-16
      &  3.1641167151_wp,  1.4177781099_wp,  4.2825156987_wp,  4.0979720404_wp, & !17-20
      &  3.4027683492_wp,  3.1330930644_wp,  3.2083129359_wp,  3.1971240284_wp, & !21-24
      &  3.0122827243_wp,  3.0215412255_wp,  2.9286697665_wp,  2.9318983659_wp, & !25-28
      &  2.9333228012_wp,  3.4211414769_wp,  3.5543149265_wp,  3.3547895521_wp, & !29-32
      &  3.8566746758_wp,  4.0522579752_wp,  3.7055903624_wp,  2.1533955559_wp, & !33-36
      &  4.8750192244_wp,  4.2251415193_wp,  3.8193395754_wp,  3.7700784196_wp, & !37-40
      &  3.8026660286_wp,  3.4791250751_wp,  3.3748738252_wp,  3.3400607232_wp, & !41-44
      &  3.3194948126_wp,  3.5185381046_wp,  3.6974620558_wp,  4.2120386946_wp, & !45-48
      &  4.2834967376_wp,  4.0408029917_wp,  4.1029792717_wp,  4.5056357496_wp, & !49-52
      &  4.1912939737_wp,  3.1889722321_wp,  5.3761906399_wp,  4.9848540155_wp, & !53-56
      &  4.1643020686_wp,  4.2242687055_wp,  4.0906998457_wp,  4.0483164017_wp, & !57-60
      &  4.0130748483_wp,  3.6618368303_wp,  3.8161213688_wp,  3.6044393411_wp, & !61-64
      &  3.7159335631_wp,  3.8610077243_wp,  3.8543967507_wp,  3.7804332520_wp, & !65-68
      &  3.6171475823_wp,  3.6614908934_wp,  3.9127576452_wp,  3.7447075110_wp, & !69-72
      &  3.7737132271_wp,  3.3371881773_wp,  3.3105897209_wp,  3.3868061092_wp, & !73-76
      &  3.4036207674_wp,  3.5310959808_wp,  3.5697281366_wp,  4.3942379403_wp, & !77-80
      &  4.6313791191_wp,  4.3952892419_wp,  4.2961541488_wp,  4.6294486870_wp, & !81-84
      &  4.5547581691_wp,  3.6325616087_wp,  5.0182872162_wp,  4.4284455579_wp, & !85-88
      &  3.7478960318_wp,  2.9868700652_wp,  3.6306659286_wp,  3.8514208797_wp, & !89-92
      &  3.4838982283_wp,  3.5107992105_wp,  3.4255592145_wp,  3.5581888894_wp, & !93-96
      &  3.3079867688_wp,  3.4972376288_wp,  3.4590842861_wp,  3.2327976004_wp, & !97-100
      &  3.4619632872_wp,  3.7360296053_wp,  3.5692246969_wp] !101-103

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
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !1-4
      &  1.0200000000_wp,  1.0000000000_wp,  1.0500000000_wp,  1.0900000000_wp, & !5-8
      &  1.1500000000_wp,  0.8700000000_wp,  1.1500000000_wp,  1.1500000000_wp, & !9-12
      &  1.0700000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !13-16
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !17-20
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !21-24
      &  1.0000000000_wp,  1.1000000000_wp,  1.1000000000_wp,  1.1000000000_wp, & !25-28
      &  1.0000000000_wp,  1.1000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !29-32
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !33-36
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !37-40
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !41-44
      &  1.1000000000_wp,  1.1000000000_wp,  1.1000000000_wp,  1.1000000000_wp, & !45-48
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !49-52
      &  1.0000000000_wp,  1.0000000000_wp,  0.8500000000_wp,  1.0000000000_wp, & !53-56
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !57-60
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !61-64
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !65-68
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !69-72
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !73-76
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.1000000000_wp, & !77-80
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !81-84
      &  1.0000000000_wp,  1.0000000000_wp,  0.8500000000_wp,  1.0000000000_wp, & !85-88
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
