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
      & get_eeqbc_kqchi, get_eeqbc_kqeta, get_eeqbc_kcnrad, get_eeqbc_cap, &
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

   !> Element-specific CN scaling of the charge widths for the EEQ_BC charges.
   interface get_eeqbc_kcnrad
      module procedure :: get_eeqbc_kcnrad_sym
      module procedure :: get_eeqbc_kcnrad_num
   end interface get_eeqbc_kcnrad

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
      &  0.8467604300_wp, -0.2183140421_wp, -0.2075642802_wp,  0.6522721704_wp, & !1-4
      &  0.9442682622_wp,  1.9804836854_wp,  2.5851212393_wp,  2.8899963829_wp, & !5-8
      &  2.5027257885_wp,  0.5948378668_wp, -0.8832496464_wp, -0.0739757977_wp, & !9-12
      &  0.0621943256_wp,  0.8438247314_wp,  1.5618039208_wp,  2.1734204952_wp, & !13-16
      &  2.4947520773_wp,  0.3027920207_wp, -1.8035382891_wp, -0.3011848531_wp, & !17-20
      & -0.5247459665_wp, -0.0522296680_wp,  0.0350320489_wp,  0.1261736863_wp, & !21-24
      &  0.0879135815_wp,  0.1566604025_wp,  0.0859505026_wp,  0.6615688577_wp, & !25-28
      &  0.3291294369_wp, -0.0128862111_wp,  0.2031753484_wp,  0.9199979851_wp, & !29-32
      &  1.4243885459_wp,  2.6521595171_wp,  1.9077722764_wp,  0.2246007577_wp, & !33-36
      & -1.1562222910_wp, -0.6934477379_wp, -0.4503077880_wp, -0.1476104597_wp, & !37-40
      & -0.1471643609_wp,  0.5015056348_wp,  0.5635437343_wp,  0.4136700018_wp, & !41-44
      &  0.3828956278_wp,  1.0306199982_wp, -0.0232888362_wp, -0.0704743109_wp, & !45-48
      &  0.1331506404_wp,  0.7757110023_wp,  0.8477433039_wp,  1.7958803155_wp, & !49-52
      &  1.7728657743_wp,  0.2535515446_wp, -1.4209559908_wp, -0.6223317850_wp, & !53-56
      & -0.9112942124_wp, -0.9912489704_wp, -1.7587751013_wp, -1.4095656646_wp, & !57-60
      & -1.6004374597_wp, -2.0746083045_wp, -1.7243963931_wp, -1.4371945362_wp, & !61-64
      & -1.4657711642_wp, -1.4075819085_wp, -1.2135148483_wp, -1.6180621446_wp, & !65-68
      & -1.5122836483_wp, -1.3295933777_wp, -1.2966952125_wp, -0.2223774675_wp, & !69-72
      &  0.2796591914_wp,  0.6579493095_wp,  0.9324788787_wp,  0.8893700791_wp, & !73-76
      &  0.7789798354_wp,  1.1716662587_wp,  0.6106974919_wp,  0.2622317365_wp, & !77-80
      &  0.2462209329_wp,  0.1130744929_wp,  0.5671442519_wp,  1.3272102619_wp, & !81-84
      &  1.6281355036_wp,  0.0873378501_wp, -1.3247454096_wp, -0.4173267854_wp, & !85-88
      & -0.7785494394_wp, -0.6099356260_wp, -0.8327423998_wp, -1.0441843860_wp, & !89-92
      & -0.5401572596_wp, -0.8144566011_wp, -0.7154181149_wp, -0.7102059966_wp, & !93-96
      & -0.5048909376_wp, -0.4963167616_wp, -0.7233431134_wp, -0.8625643907_wp, & !97-100
      & -0.6411108832_wp, -0.8832896019_wp, -1.0659105468_wp] !101-103


   !> Element-specific chemical hardnesses for the EEQ_BC charges.
   real(wp), parameter :: eeqbc_eta(max_elem) = [&
      &  4.5883961733_wp,  1.1957515220_wp, -0.2947755074_wp, -0.1757316605_wp, & !1-4
      & -0.7242416016_wp, -0.8846215047_wp, -1.2478764437_wp, -0.8814120281_wp, & !5-8
      & -2.8812974754_wp, -0.7235501142_wp, -2.9340352578_wp, -2.2205741355_wp, & !9-12
      & -2.5267711878_wp, -1.5191866025_wp, -1.7594603150_wp, -1.0444104706_wp, & !13-16
      & -0.6593014736_wp, -1.8072342920_wp,  0.9319470630_wp, -0.4716460928_wp, & !17-20
      & -0.4810207743_wp, -0.5998907116_wp, -1.0259305973_wp, -0.5232068685_wp, & !21-24
      & -0.8438670173_wp, -1.4136461445_wp, -1.4302672869_wp, -0.3590539544_wp, & !25-28
      & -0.2294329581_wp, -0.2640773231_wp, -0.3079907817_wp, -1.5584979933_wp, & !29-32
      & -1.3523584209_wp, -0.4009220841_wp,  0.8019234081_wp,  1.2644457366_wp, & !33-36
      &  0.4973132401_wp, -1.9189725404_wp, -1.2158180734_wp, -0.5850266122_wp, & !37-40
      & -0.1577964029_wp,  0.0695444240_wp,  0.1222721193_wp,  0.1936223827_wp, & !41-44
      &  0.2215489239_wp,  0.0883632091_wp,  0.5002471655_wp, -1.2704295500_wp, & !45-48
      & -1.1831319857_wp, -1.6036372808_wp, -1.4774758423_wp, -0.3040420067_wp, & !49-52
      &  0.9836010535_wp,  1.2274190712_wp,  0.4713569859_wp, -1.2186383438_wp, & !53-56
      & -1.1159875541_wp, -2.3297550997_wp, -1.5576519828_wp, -1.5338533611_wp, & !57-60
      & -1.6850521824_wp, -1.7849259848_wp, -1.6940014126_wp, -1.3584868519_wp, & !61-64
      & -1.9445485275_wp, -2.4592398168_wp, -1.8480509365_wp, -1.5205053566_wp, & !65-68
      & -0.4655382999_wp, -0.7559411426_wp, -0.5314313657_wp, -0.6771778094_wp, & !69-72
      & -0.4595762663_wp, -0.1546279733_wp,  0.1400396898_wp,  0.3053897273_wp, & !73-76
      &  0.2952195420_wp,  0.2218883373_wp,  0.1605445237_wp,  0.1077340830_wp, & !77-80
      & -0.2375633028_wp, -0.3386633415_wp, -1.2574138276_wp, -0.9345365231_wp, & !81-84
      & -0.6371937508_wp,  0.3119375461_wp, -0.1411132568_wp, -3.3424645158_wp, & !85-88
      & -0.8506132555_wp, -0.8565428723_wp, -0.7712832848_wp, -0.5526607160_wp, & !89-92
      & -0.5472725142_wp, -0.5673325711_wp, -0.5244494931_wp, -0.4080818111_wp, & !93-96
      & -0.6009383615_wp, -1.0596382731_wp, -0.9398333072_wp, -0.7472850212_wp, & !97-100
      & -1.0869214525_wp, -1.5353612896_wp, -2.2894932658_wp] !101-103

   !> Element-specific charge widths for the EEQ_BC charges.
   real(wp), parameter :: eeqbc_rad(max_elem) = [&
      &  4.0000000000_wp,  0.0355969694_wp,  0.2650154420_wp,  0.5389107261_wp, & !1-4
      &  0.3282695634_wp,  0.2651875667_wp,  0.2554623021_wp,  0.3013517866_wp, & !5-8
      &  0.1637327020_wp,  0.0387183510_wp,  0.1525511560_wp,  0.2155124716_wp, & !9-12
      &  0.1692980040_wp,  0.2546553636_wp,  0.2250243725_wp,  0.2682171535_wp, & !13-16
      &  0.2468404216_wp,  0.1529624346_wp,  0.6885263237_wp,  0.4768491463_wp, & !17-20
      &  0.4078987025_wp,  0.4306915507_wp,  0.3421409676_wp,  0.4995155317_wp, & !21-24
      &  0.4856480110_wp,  0.3717869315_wp,  0.3037240219_wp,  0.5974648434_wp, & !25-28
      &  0.4342799350_wp,  0.3980074086_wp,  0.3421199189_wp,  0.2879133984_wp, & !29-32
      &  0.2996040742_wp,  0.2923341252_wp,  0.5785411499_wp,  0.5997355731_wp, & !33-36
      &  0.4693960799_wp,  0.2740103366_wp,  0.3607234702_wp,  0.4212613141_wp, & !37-40
      &  0.4684912570_wp,  0.7362036201_wp,  0.9950778239_wp,  0.7472159231_wp, & !41-44
      &  0.7873394714_wp,  0.8313984448_wp,  0.4442011936_wp,  0.2641990747_wp, & !45-48
      &  0.3039104080_wp,  0.2891944683_wp,  0.2924757198_wp,  0.4412748805_wp, & !49-52
      &  0.6920797463_wp,  0.7728476203_wp,  0.5170524845_wp,  0.4105460865_wp, & !53-56
      &  0.3741823773_wp,  0.2269034230_wp,  0.2427137472_wp,  0.2532377090_wp, & !57-60
      &  0.2554764426_wp,  0.2480417102_wp,  0.2370393810_wp,  0.2728174166_wp, & !61-64
      &  0.2270205599_wp,  0.2027271571_wp,  0.2210537791_wp,  0.2346159975_wp, & !65-68
      &  0.3174302311_wp,  0.2701589120_wp,  0.3262664264_wp,  0.4177075626_wp, & !69-72
      &  0.5303517152_wp,  0.7353995137_wp,  1.0181499685_wp,  0.9916888957_wp, & !73-76
      &  0.9204759171_wp,  0.6599501949_wp,  0.7116202062_wp,  0.4727757004_wp, & !77-80
      &  0.4137738690_wp,  0.5004528005_wp,  0.3432877772_wp,  0.3523498261_wp, & !81-84
      &  0.3662499879_wp,  0.4212391629_wp,  0.2392397490_wp,  0.1864832533_wp, & !85-88
      &  0.4664056499_wp,  0.5143416354_wp,  0.4921261258_wp,  0.5206678560_wp, & !89-92
      &  0.5796701081_wp,  0.5083748828_wp,  0.5103017400_wp,  0.5269057811_wp, & !93-96
      &  0.4779200223_wp,  0.3674815329_wp,  0.3428177247_wp,  0.4392766993_wp, & !97-100
      &  0.3364046219_wp,  0.2608548259_wp,  0.2227971082_wp] !101-103

   !> Element-specific CN scaling of the electronegativity for the EEQ_BC charges.
   real(wp), parameter :: eeqbc_kcnchi(max_elem) = [&
      &  5.5328922696_wp, 12.9737820266_wp, -1.6592106860_wp, -0.3154196961_wp, & !1-4
      &  0.4441249481_wp,  1.5896959149_wp,  1.6856817044_wp,  2.4835888029_wp, & !5-8
      &  3.7064657438_wp, 18.1918276402_wp, -1.5817634740_wp, -0.7114945163_wp, & !9-12
      &  3.0113116005_wp,  1.0430409506_wp,  1.4052802324_wp,  2.1644738944_wp, & !13-16
      &  6.0232628785_wp, 27.8714382497_wp, -4.7499672972_wp, -2.2386395673_wp, & !17-20
      & -1.4252904509_wp, -1.9542868863_wp, -1.5285798337_wp, -2.2844491884_wp, & !21-24
      & -2.4931928871_wp, -2.1928513769_wp, -1.8066244309_wp, -1.6257495871_wp, & !25-28
      & -0.4640973167_wp, -0.3441517545_wp,  0.0793393574_wp,  0.0676529583_wp, & !29-32
      &  1.0499008514_wp,  3.5218993653_wp,  3.8194162606_wp,  7.0195831367_wp, & !33-36
      & -2.4670763641_wp, -1.5229101143_wp, -1.2573177520_wp, -1.5811870823_wp, & !37-40
      & -1.2957183160_wp, -0.5867834716_wp, -1.4501119192_wp, -1.1263505554_wp, & !41-44
      & -0.9577820087_wp, -0.0560379401_wp,  0.1385839494_wp,  0.8184372572_wp, & !45-48
      &  0.3593414234_wp,  0.7072574068_wp,  0.5700824894_wp,  2.1968779757_wp, & !49-52
      &  3.8579669467_wp,  5.0467400400_wp, -5.9921627529_wp, -2.6864064013_wp, & !53-56
      & -1.9441641042_wp, -1.3844301600_wp, -2.7597617890_wp, -2.1393765148_wp, & !57-60
      & -2.4624820159_wp, -2.4578802325_wp, -1.9281082277_wp, -1.8867412974_wp, & !61-64
      & -2.0320133577_wp, -1.8933320845_wp, -1.6235701523_wp, -1.6660649802_wp, & !65-68
      & -1.2972155967_wp, -1.4577593511_wp, -1.1549925150_wp, -2.1993827525_wp, & !69-72
      & -1.2084619352_wp, -0.3372745183_wp, -0.8104435845_wp, -0.7919301473_wp, & !73-76
      & -0.7662020648_wp, -0.0827687602_wp,  0.7765541362_wp,  0.0846993657_wp, & !77-80
      &  0.0658165860_wp,  0.0664742859_wp,  0.0717231616_wp,  1.4933045369_wp, & !81-84
      &  2.9655886641_wp,  6.8243798237_wp, -5.1925813067_wp, -2.8796122621_wp, & !85-88
      & -5.3726978306_wp, -9.0166913689_wp, -3.6051413102_wp, -3.6532963392_wp, & !89-92
      & -3.5589368234_wp, -2.6139340025_wp, -3.1522330210_wp, -2.6905611959_wp, & !93-96
      & -1.9788157098_wp, -1.8145174473_wp, -1.6079639366_wp, -2.8677464047_wp, & !97-100
      & -1.5325075455_wp, -0.9706879719_wp, -0.9870924795_wp] !101-103

   !> Element-specific local q scaling of the electronegativity for the EEQ_BC charges.
   real(wp), parameter :: eeqbc_kqchi(max_elem) = [&
      &  3.5194720929_wp,  0.1470383338_wp, 11.3046142065_wp,  7.5532876529_wp, & !1-4
      &  6.9405502860_wp,  5.3352972531_wp,  5.6154161978_wp,  5.3451107586_wp, & !5-8
      &  4.9238939530_wp,  0.3010529968_wp, 11.4114585448_wp,  8.7529247902_wp, & !9-12
      &  8.4595779428_wp,  8.3601655718_wp,  6.8560006460_wp,  6.9560862453_wp, & !13-16
      &  6.5872746910_wp,  1.0284627215_wp, 11.2332729749_wp,  8.3518727777_wp, & !17-20
      &  8.2180840237_wp,  8.2081092109_wp,  5.8720389297_wp,  6.2939978919_wp, & !21-24
      &  5.2263267105_wp,  6.3757018080_wp,  7.4588049420_wp,  7.5537205523_wp, & !25-28
      &  7.0423424647_wp,  7.7917037233_wp,  8.9982887212_wp,  7.9158554904_wp, & !29-32
      &  7.4426719387_wp,  7.9251233756_wp,  7.4010883235_wp,  2.0437940547_wp, & !33-36
      &  9.7654126363_wp,  8.6081160042_wp,  8.1112159520_wp,  8.1587876633_wp, & !37-40
      &  7.0493777770_wp,  5.9157675086_wp,  5.9803404959_wp,  6.3483633640_wp, & !41-44
      &  6.5377419233_wp,  6.9235703046_wp,  7.1002092982_wp,  7.6951020528_wp, & !45-48
      &  9.0456773925_wp,  8.1988149793_wp,  8.7972736262_wp,  7.2599240745_wp, & !49-52
      &  8.3749337039_wp,  3.0132713616_wp, 11.5073160210_wp,  9.6129336117_wp, & !53-56
      &  8.1161664515_wp,  7.2799080442_wp,  5.8097665286_wp,  7.0148831513_wp, & !57-60
      &  7.5726441117_wp,  6.5689836152_wp,  7.2364344076_wp,  7.3499748371_wp, & !61-64
      &  6.9537454298_wp,  7.6680778361_wp,  7.8915854710_wp,  7.3663222748_wp, & !65-68
      &  8.3295557281_wp,  9.1537837748_wp,  8.3384155195_wp,  7.7962587778_wp, & !69-72
      &  7.6219560496_wp,  6.5513950643_wp,  5.8552610330_wp,  6.5936071685_wp, & !73-76
      &  7.1151593969_wp,  7.5468614625_wp,  7.4826342310_wp,  7.4331660695_wp, & !77-80
      &  9.3548144042_wp,  8.6417563154_wp,  7.7980785441_wp,  8.3455968876_wp, & !81-84
      &  8.0533780237_wp,  2.7174467656_wp,  8.6880898876_wp,  9.7678340910_wp, & !85-88
      &  7.1580873338_wp,  6.1870983937_wp,  6.7235256040_wp,  6.3349783291_wp, & !89-92
      &  7.6603424102_wp,  7.5825328428_wp,  7.5569480689_wp,  7.0155860477_wp, & !93-96
      &  7.6040130324_wp,  9.3179608505_wp,  8.0442848320_wp,  7.5793455390_wp, & !97-100
      &  7.8843649728_wp,  7.9778579082_wp,  7.4274418789_wp] !101-103

   !> Element-specific local q scaling of the chemical hardness for the EEQ_BC charges
   real(wp), parameter :: eeqbc_kqeta(max_elem) = 0.5_wp * [&
      &  7.5990940066_wp, -1.7712337517_wp,  6.1461685618_wp,  2.8690613059_wp, & !1-4
      &  1.6643432851_wp,  0.5232151629_wp,  0.4094442242_wp,  0.5445440043_wp, & !5-8
      &  0.5602960536_wp, -3.6652419435_wp,  7.7726999313_wp,  1.7981461482_wp, & !9-12
      &  1.6479386985_wp,  1.5352294048_wp,  1.4476890514_wp,  0.3873967330_wp, & !13-16
      &  0.3447149070_wp, -3.7066804134_wp, 10.0299413574_wp,  0.8409780941_wp, & !17-20
      &  0.6998790974_wp,  1.2389616714_wp,  0.2201436209_wp,  0.1049052716_wp, & !21-24
      &  0.1701168098_wp,  0.7729753145_wp,  0.9625653597_wp,  0.8681267002_wp, & !25-28
      &  1.0664746178_wp,  1.1490788354_wp,  0.8323147992_wp,  0.1259877101_wp, & !29-32
      &  1.8338867331_wp,  1.3409108611_wp,  1.7417458501_wp, -3.0797688943_wp, & !33-36
      &  4.5253732961_wp,  1.6957541047_wp,  0.7089513921_wp,  0.6159555835_wp, & !37-40
      &  0.1997696824_wp,  0.2894247478_wp,  0.6439800012_wp,  0.6100768786_wp, & !41-44
      &  0.5265676353_wp,  0.2620956605_wp,  0.4528880198_wp,  0.9150343377_wp, & !45-48
      &  1.1181884701_wp,  1.1416131571_wp,  2.5598947488_wp,  0.7705458485_wp, & !49-52
      &  2.0530695641_wp, -3.3633902941_wp,  6.8090216039_wp,  1.8699226003_wp, & !53-56
      &  0.4451894480_wp,  0.6210853454_wp,  0.0336308783_wp,  0.2286470656_wp, & !57-60
      &  0.4072821694_wp,  0.5218779999_wp,  0.4608098284_wp,  0.3882996289_wp, & !61-64
      &  0.3072542765_wp,  0.3874470962_wp,  0.4505921870_wp,  0.3947682979_wp, & !65-68
      &  0.9675164017_wp,  0.7942670664_wp,  0.8496733461_wp,  0.1959035765_wp, & !69-72
      &  0.8754101873_wp,  0.6473023504_wp,  0.2780192518_wp,  0.3030047206_wp, & !73-76
      &  1.1869873530_wp,  1.1081891092_wp,  1.0776929512_wp,  1.2264041160_wp, & !77-80
      &  1.2040659543_wp,  1.4408705128_wp,  1.6466651020_wp,  1.0202422264_wp, & !81-84
      &  1.0468694472_wp, -3.1816369086_wp, -3.8817659859_wp,  3.4291745704_wp, & !85-88
      & -0.1472221610_wp, -0.2871652305_wp, -0.1594027016_wp, -0.3466228137_wp, & !89-92
      &  0.1939883436_wp,  0.9013047150_wp,  0.6057392725_wp,  0.1227569526_wp, & !93-96
      &  0.2399047925_wp,  1.6413012417_wp,  0.7875024434_wp,  0.3350893161_wp, & !97-100
      &  0.6812462134_wp,  0.3920559029_wp,  0.1168443747_wp] !101-103

   !> Element-specific CN scaling of the charge widths for the EEQ_BC charges.
   real(wp), parameter :: eeqbc_kcnrad(max_elem) = [&
      &  0.5988208600_wp,  0.3321474006_wp,  0.4316952733_wp,  0.4104987934_wp, & !1-4
      &  0.3690236763_wp,  0.4270347634_wp,  0.5320217767_wp,  0.5725076800_wp, & !5-8
      &  0.3940113797_wp,  0.3736658257_wp,  0.2864468058_wp,  0.2723805376_wp, & !9-12
      &  0.2491105505_wp,  0.2755612552_wp,  0.2910633351_wp,  0.3475200063_wp, & !13-16
      &  0.3098214060_wp,  0.2864468058_wp,  0.3609976027_wp,  0.3840352955_wp, & !17-20
      &  0.3085826265_wp,  0.3020537093_wp,  0.4010007764_wp,  0.3897192745_wp, & !21-24
      &  0.3356093778_wp,  0.3561348788_wp,  0.3289902150_wp,  0.3332060568_wp, & !25-28
      &  0.3349923518_wp,  0.2491105505_wp,  0.2491105505_wp,  0.2853356568_wp, & !29-32
      &  0.2621597242_wp,  0.3157754405_wp,  0.3573870228_wp,  0.2906289756_wp, & !33-36
      &  0.3425107254_wp,  0.3319662046_wp,  0.4043812044_wp,  0.2634137659_wp, & !37-40
      &  0.2770233690_wp,  0.2978081418_wp,  0.2491105505_wp,  0.2544302365_wp, & !41-44
      &  0.3519334025_wp,  0.3563832964_wp,  0.2906289756_wp,  0.2491105505_wp, & !45-48
      &  0.2751304809_wp,  0.2706509751_wp,  0.2566203511_wp,  0.3216109243_wp, & !49-52
      &  0.3260885700_wp,  0.2906289756_wp,  0.2974124210_wp,  0.3885374374_wp, & !53-56
      &  0.2717531659_wp,  0.3412104084_wp, -0.2446020708_wp,  0.2716777886_wp, & !57-60
      &  0.2996438556_wp,  0.2782260681_wp,  0.3031591013_wp,  0.2914395673_wp, & !61-64
      &  0.2491105505_wp,  0.3068333755_wp,  0.2530102617_wp,  0.2783134073_wp, & !65-68
      &  0.2809577441_wp,  0.2491105505_wp,  0.2491105505_wp,  0.3423325847_wp, & !69-72
      &  0.3652718437_wp,  0.2575213346_wp,  0.2491105505_wp,  0.2491105505_wp, & !73-76
      &  0.3117867169_wp,  0.2615212840_wp,  0.2433103566_wp,  0.2491105505_wp, & !77-80
      &  0.3267638312_wp,  0.2455258336_wp,  0.2602384046_wp,  0.2658671290_wp, & !81-84
      &  0.2444132535_wp,  0.2455258336_wp,  0.3059735459_wp,  0.2455258336_wp, & !85-88
      &  0.2721190118_wp,  0.3097788034_wp,  0.2530183420_wp,  0.3171856474_wp, & !89-92
      &  0.2928452400_wp,  0.3509598360_wp,  0.3289437683_wp,  0.3243239628_wp, & !93-96
      &  0.3506262769_wp,  0.3049422816_wp,  0.2491105505_wp,  0.2785530744_wp, & !97-100
      &  0.3191189832_wp,  0.3516748198_wp,  0.2918605964_wp] !101-103

   !> Element-specific bond capacitance for the EEQ_BC charges.
   real(wp), parameter :: eeqbc_cap(max_elem) = [&
      &  0.3428467448_wp,  0.1106892417_wp,  0.3514036826_wp,  0.7569234600_wp, & !1-4
      &  0.7270507604_wp,  1.1490434838_wp,  0.8694873780_wp,  0.8380832427_wp, & !5-8
      &  1.1244986101_wp,  0.4919068115_wp,  0.3159074795_wp,  0.4776472442_wp, & !9-12
      &  0.2221030280_wp,  0.5105140282_wp,  0.7826588790_wp,  0.7390711438_wp, & !13-16
      &  0.7974090328_wp,  0.2249646212_wp,  0.2578475367_wp,  0.7390346797_wp, & !17-20
      &  0.5591334361_wp,  0.5465652686_wp,  0.6245310548_wp,  1.4259897704_wp, & !21-24
      &  1.2214664278_wp,  1.5636510391_wp,  1.4618416034_wp,  1.2025691607_wp, & !25-28
      &  1.0987721637_wp,  0.6909015748_wp,  0.5459114123_wp,  0.6295036380_wp, & !29-32
      &  1.0027509834_wp,  0.5985583574_wp,  1.5487838126_wp,  0.2507548457_wp, & !33-36
      &  0.6415065092_wp,  0.5309302633_wp,  0.7044774615_wp,  0.3809101340_wp, & !37-40
      &  0.7395660266_wp,  1.3806539570_wp,  1.3796773107_wp,  1.6171291196_wp, & !41-44
      &  1.9236192491_wp,  2.3300214417_wp,  0.6466145903_wp,  0.5084074102_wp, & !45-48
      &  0.6661541273_wp,  0.8948242500_wp,  0.9027397906_wp,  0.7044283953_wp, & !49-52
      &  1.4716963736_wp,  0.3244580156_wp,  0.8963931782_wp,  0.6553611168_wp, & !53-56
      &  0.6640754386_wp,  0.1929343121_wp,  0.2419963193_wp,  0.3056746183_wp, & !57-60
      &  0.2045621079_wp,  0.6162340144_wp,  0.2964476457_wp,  0.3282269293_wp, & !61-64
      &  0.4431835444_wp,  0.3799918471_wp,  0.4070974521_wp,  0.3616619272_wp, & !65-68
      &  0.3810084239_wp,  0.3361447153_wp,  0.2885834022_wp,  0.6530717463_wp, & !69-72
      &  0.6453994519_wp,  1.4288339489_wp,  1.4831346363_wp,  1.2919121889_wp, & !73-76
      &  1.5535370198_wp,  2.6046063268_wp,  2.0913235323_wp,  1.0601341622_wp, & !77-80
      &  0.7255818048_wp,  1.1334579037_wp,  1.4645026198_wp,  0.9794476310_wp, & !81-84
      &  1.1741195488_wp,  0.2719384127_wp,  0.2621496150_wp,  0.7013594217_wp, & !85-88
      &  0.5078479860_wp,  0.3916895939_wp,  0.3312179340_wp,  0.4097939927_wp, & !89-92
      &  0.4927298461_wp,  0.4132140952_wp,  0.5534174462_wp,  0.5716238505_wp, & !93-96
      &  0.6090670120_wp,  0.4414188794_wp,  0.5416479318_wp,  0.6434739754_wp, & !97-100
      &  0.3827783674_wp,  0.3287073050_wp,  0.5096472311_wp] !101-103

   !> Element-specific covalent radii for the CN for the EEQ_BC charges.
   real(wp), parameter :: eeqbc_cov_radii(max_elem) = 0.5_wp*[&
      &  0.5930850061_wp,  0.6828783524_wp,  2.5060627321_wp,  2.1530174721_wp, & !1-4
      &  2.0808370973_wp,  2.1547612587_wp,  2.3825166136_wp,  2.3845548778_wp, & !5-8
      &  2.2785689025_wp,  1.5895245890_wp,  3.2262348700_wp,  3.0085723165_wp, & !9-12
      &  2.4255063979_wp,  2.8239473425_wp,  3.1186237627_wp,  3.3331763174_wp, & !13-16
      &  3.1216907679_wp,  1.5732371114_wp,  3.2900208529_wp,  3.4202206471_wp, & !17-20
      &  2.9200724722_wp,  3.0066744258_wp,  2.8630113712_wp,  2.8148129713_wp, & !21-24
      &  2.4967937635_wp,  3.0257979292_wp,  3.1144574229_wp,  2.8415686668_wp, & !25-28
      &  2.3076713228_wp,  2.7929299468_wp,  2.9384555000_wp,  3.1211972610_wp, & !29-32
      &  3.3019239452_wp,  3.5482037468_wp,  3.5268705061_wp,  2.3017605635_wp, & !33-36
      &  3.6453427844_wp,  3.6627725418_wp,  3.1732780680_wp,  3.3336121647_wp, & !37-40
      &  3.1865126417_wp,  3.1583027293_wp,  2.9482110819_wp,  3.1853598435_wp, & !41-44
      &  3.2430942807_wp,  3.1951785079_wp,  2.9793621882_wp,  3.1543875559_wp, & !45-48
      &  3.1850367978_wp,  3.5529181164_wp,  3.4500895536_wp,  3.8418786490_wp, & !49-52
      &  3.8953057563_wp,  2.8735372659_wp,  3.7155011702_wp,  3.9719089931_wp, & !53-56
      &  3.1569703924_wp,  3.4936913809_wp,  3.3335208117_wp,  3.9243845028_wp, & !57-60
      &  3.8732096780_wp,  3.1689369113_wp,  3.4987436252_wp,  3.5190884141_wp, & !61-64
      &  3.6663276643_wp,  3.7941855096_wp,  3.6947440641_wp,  3.7383619093_wp, & !65-68
      &  3.7065668579_wp,  3.5989802683_wp,  3.5113445936_wp,  3.3622392888_wp, & !69-72
      &  3.0067203462_wp,  3.0264024510_wp,  2.8716706964_wp,  3.2964048554_wp, & !73-76
      &  3.3131562414_wp,  3.4380505276_wp,  2.9801688586_wp,  3.6006573663_wp, & !77-80
      &  3.4878494275_wp,  4.0080167544_wp,  3.4867375894_wp,  3.8412537998_wp, & !81-84
      &  4.0986699510_wp,  3.0264945260_wp,  4.1292282410_wp,  3.7772409489_wp, & !85-88
      &  1.8907313528_wp,  1.0395376590_wp,  1.9323562902_wp,  3.0573080720_wp, & !89-92
      &  2.7107637944_wp,  3.3756041761_wp,  3.2094752255_wp,  3.2207515613_wp, & !93-96
      &  2.7188494129_wp,  2.7941335957_wp,  3.2020860461_wp,  1.5623494196_wp, & !97-100
      &  2.8208057454_wp,  2.8130190589_wp,  2.7482546063_wp] !101-103

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
      &  1.0300000000_wp,  1.0000000000_wp,  1.0500000000_wp,  1.0200000000_wp, & !5-8
      &  1.1800000000_wp,  0.8500000000_wp,  1.0500000000_wp,  1.1500000000_wp, & !9-12
      &  1.1000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !13-16
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !17-20
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !21-24
      &  1.0000000000_wp,  1.1000000000_wp,  1.1000000000_wp,  1.1000000000_wp, & !25-28
      &  1.0000000000_wp,  1.1500000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !29-32
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !33-36
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !37-40
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !41-44
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.1500000000_wp, & !45-48
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !49-52
      &  1.0000000000_wp,  1.0000000000_wp,  0.5000000000_wp,  1.0000000000_wp, & !53-56
      &  1.2000000000_wp,  1.2000000000_wp,  1.2000000000_wp,  1.2000000000_wp, & !57-60
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

!> Get CN scaling of the charge width for species with a given symbol
elemental function get_eeqbc_kcnrad_sym(symbol) result(kcnrad)

   !> Element symbol
   character(len=*), intent(in) :: symbol

   !> CN scaling of charge width
   real(wp) :: kcnrad

   kcnrad = get_eeqbc_kcnrad(to_number(symbol))

end function get_eeqbc_kcnrad_sym

!> Get CN scaling of the charge width for species with a given atomic number
elemental function get_eeqbc_kcnrad_num(number) result(kcnrad)

   !> Atomic number
   integer, intent(in) :: number

   !> CN scaling of charge width
   real(wp) :: kcnrad

   if (number > 0 .and. number <= size(eeqbc_kcnrad, dim=1)) then
      kcnrad = eeqbc_kcnrad(number)
   else
      kcnrad = -1.0_wp
   end if

end function get_eeqbc_kcnrad_num

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
