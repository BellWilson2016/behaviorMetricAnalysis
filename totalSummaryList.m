function tSL = totalSummaryList()

	% Up to date through 4036 plot #133
	% Removed 3641, 3651, 3666, 3669, 3671, 3676, 3677, 3679, 3685, 3745 due to
	%	NorpA ; H134R parent stock contamination. (7/1 - 8/15)
	% Remove 3704 due to air error, 3694 due to laser arm error
	%
	tSL =      {{'Or7a II',3517,3518,3567,3571},... 
				{'Or7a III',3484,3487}, ...
				{'Or10a',3493,3494},...
				{'Gr21a II Suh',3481,3483,3537},...
				{'Gr21a III BD1',3507,3512},...
				{'Gr21a III BD2',3496}...
				{'Or22a',3482,3485,3538,3540,3721,3887,3889},...
				{'Or42a',3490,3547,3550,3742,3905,3947},...
				{'Or42b',3474,3475,3530,3533,3698,3953},...
				{'Or56a II',3516,3519,3823,3886,3890},...
				{'Or56a III',3491,3492,3633,3827,3883},...
				{'Gr63a II',3508,3566,3569,3753,3879},...
				{'Gr63a III',3498,3499,3822,3877,3897,3900},...
				{'Ir64a III',3513,3615,3618,3829,3832,3909},...
				{'Or67b',3502,3506,3621,3624,3857,3860},...
				{'Or67d',3488,3539,3542,3636,3862},...
				{'Or82a',3495,3497,3543}...
				{'Or83b',3479,3480,3532,3534,3798,3804,3808,3814,3816},...
				{'Or85a II',3522,3525,3622,3846,3916,3917,3919},...
				{'Or85a III',3528,3531,3731,3912,3914},...
				{'Or85b II',3520,3524},...
				{'Or85b III',3514,3515,3584},...
				{'Or92a',3486,3489,3535,3536,3836,3838},...
				{'Ctrl',3476,3477,3478,3500,3548,3637,3747,3888},...
				{'Or42b ; Or22a',3541,3544,3602,3754},...
				{'Bl ; Or22a / Gr63a',3545,3657,3784,3815,3817},...
				{'Or67d ; Or22a', 3546,3552,3713,3716},...
				{'Or42b ; Or92a', 3549,3551,3606,3672},...
				{'Or42a ; Or22a', 3553,3555,3695,3774},...
				{'Bl ; Or92a / Gr63a',3554,3598,3824,3826},...
				{'Or22a / Or92a', 3556,3558,3737,3761,3763},...
				{'Or56a ; Gr63a', 3557, 3562,3614,3720,3722},...
				{'Or56a ; Or22a', 3559, 3560, 3725, 3793},...
				{'Or85a ; Or92a',3561,3564,3696,3792},...
				{'Or85a ; Or22a',3563,3565,3599,3601},...
				{'Or42b ; Gr63a',3568,3570,3577,3782,3787,3789},...
				{'Or85a ; Gr63a',3572,3575,3692,3781,3791,3797},... 
				{'Or22a / Or82a',3573},...
				{'Or42a ; Or92a',3574,3578,3758,3834,3842},...
				{'Or42a ; Gr63a',3576,3579,3587,3623,3627},...
				{'Or22a / Or67b',3580,3583,3585,3777,3778},...
				{'Or67d ; Or92a',3581,3582,3767,3770},... 
				{'Or42b ; Ir64a',3586,3588,3594,3768},...
				{'Or92a / Or67b',3589,3593,3783,3786},...
				{'Or42b ; Or56a',3590,3592,3600,3790,3794,3796},...
				{'Or67d ; Gr63a',3595,3673,3675,3678},...
				{'Ir64a / Gr63a',3604,3616,3619,3620,3876},...
				{'Ir64a / Or92a',3608,3733,3734,3762},...
				{'Or85a ; Or56a',3612,3613,3719,3723},...
				{'42b ; 92a W-', 3610},...
				{'Or42a ; Or56a ', 3617,3691,3693,3780,3810},...
				{'H/H ; Gr63a/Gr63a',3625,3630,3799,3803},...
				{'Ir40a ; TM2',3626,3629,3833,3844,3847,3932,3936},...
				{'Or56a / Or67b', 3628,3631,3741,3743},...
				{'Or42b ; Or67b',3632,3635,3687,3690},...
				{'Or42a ; Or67b',3638,3674,3771,3773,3775},...
				{'67d HighFreq No wind',3639},... 
				{'Ir64a / Or67b', 3640,3735,3740,3785,3788},...
				{'Or67d ; Or56a',3644,3649,3706,3707},...
				{'Or85a ; Or67b',3648,3650,3712,3717},...
				{'Or67d ; Or67b',3652,3653,3703,3746},...
				{'Or56a ; Or92a',3654,3656,3800,3828,3830},...
				{'Gr63a / Or67b',3655,3658,3699,3701},...
				{'Or42a ; Ir64a',3667,3727,3748,3755,3759},...
				{'Or67d ; Ir64a',3686,3689,3805,3807,3811},...
				{'Bl ; Or56a / Ir64a',3688,3728,3729,3750},...
				{'Bl ; Or22a / Ir64a',3700,3709,3739,3801},...
				{'Ir40a ; Gr63a / TM2',3708,3710,3921,3923,3927},...
				{'Or42b ; Or85a',3711,3714,3760,3764},...
				{'Or42b III',3715,3718,3849,3855,3929,3934},...
				{'Or42a III',3724,3726,3730,3851,3854},...
				{'Or42a III Low Power',3732},...
				{'Gr63a III / Or42b III',3738,3820,3825,3837,3843},...
				{'Or67d ; Or85a',3744,3869,3872,3875},...
				{'Or83b ; Gr63a',3749,3751,3752,3812,3819},...
				{'Or22a / Or42b III',3756,3757,3831},...
				{'Or22a / Or42a III',3765,3769,3839},...
				{'Or67b / Or42a III',3766,3772},...
				{'H/H ; Or67b/Or67b',3776,3779},...
				{'H/H ; Or22a/Or22a',3806,3809},...
				{'Or83b - no wind',3802,3813,3818,3954},...
				{'Gr63a III / Or42a III',3835,3840,3930,3933},...
				{'-------'},...
				{'H/H ; Or56a/Or56a',3841},...
				{'Or85a ; Ir64a',3845,3961,3962,3966},...
				{'Or56a ; Or67b',3848,3852,3920,3922,3926},...
				{'Or56a / Or42a III',3850,3853,3880,3952},...
				{'Or92a / Or42b III',3856,3858,3991},...
				{'Or42a ; Or85a',3859,3866,3885,3943},...
				{'Gr63a ; Or92a / TM2',3861,3882,3894,3902,3908},...
				{'Or67b / Or42b III',3863,3865,3975,3979},...
				{'Or56a / Or42b III',3864,3867,3949,3974},...
				{'Or42a ; Or42b',3868,3870,3904,3985},...
				{'Gr63a ; Or56a / TM2',3871,3874,3895,3960},...
				{'Or67d ; Or42b',3878,3881,3913,3918},... 
				{'Or22a / Or85a',3884,3892,3898,3907},...
				{'Gr63a ; Or67b / TM2',3891,3896,3911,3950},...
				{'Gr63a ; Ir64a / TM2',3893,3899,3906,3910},...
				{'Or56a / Or85a',3901},...
				{'Gr63a ; Or42b / TM2',3903},...
				{'Or92a / Or85a', 3915},...
				{'Ir40a ; Or92a / TM2', 3924,3986,4001},...
				{'Ctrl - no wind',3925,3928,3993},...
				{'Ir40a ; Ir64a / TM2', 3931,3935,3937,3959},...
				{'Or42b III - LPf',3938,3939,3941,4003},...
				{'Or22a / Or92a - LPf',3940,3942,3944,3946},...
				{'Or42a III - LPf',3945,3978,3980,3982,3998},...
				{'Or92a / Or42b III - LPf',3948,3951,3957,3990,4005},...
				{'Or22a / Or42b III - LPf',3955,3956,3958,3968},...
				{'Or22a - LPf',3963,3987},...
				{'Or67d ; Or42a III',3964,3965,3967,3976},...
				{'Or92a - LPf',3970,3971,3992,3994,3997},...
				{'Gr63a III - BGO ACV',3981,3983,3984,3988,4002},...
				{'Gr63a III - Unstarved', 3989,4006,4010,4014},...
				{'Ctrl - one way wind',3996,3997},...
				{'Or42b III - one way wind',3999,4000},...
				{'Or83b - OLfifthsec',3684,3705},...
				{'Or42b II - OLfifthsec',3697},...
				{'Or42b II - OLhalfsec',3702},...
				{'Or83b - OLonesec',3105},...
				{'Gr63a III - LPf',4004,4015,4025,4030},...
				{'Ir40a II - BGO',4007,4008,4011,4013,4018},...
				{'Or83b - one way wind',4009,4012,4028,4029},...
				{'Ctrl - BGO',4016,4017,4019},...
				{'Gr63a II - BGO',4020,4021},...
				{'Ctrl - BGO LP',4022},...
				{'Ir40a II - BGO LP',4023,4031},...
				{'Gr63a II - LP',4024,4026},...
				{'Ir40a ; Or42b / TM2',4027},...
				{'Or83b - OL One Sec',4033},...
				{'Or83b - one way wind LP',4035},...
				{'Gr63a III - BGO LP',4034},...
				{'Or83b - OL Five Sec',4036},...
				};% Or42b II - OLhalfsec # 119 % Or83b - OLonesec #120 % Gr63a III - LPf #121
				% Ir40a II - BGO #122 
				% Or83b - one way wind #123
				% OR83b - one way wind LP #131
				% OR83b - OL 5 Sec #133
