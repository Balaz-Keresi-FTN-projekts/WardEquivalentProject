Ybusdom = [...
-50   0   0   40   0   0   10   0   0   0   0   0   0   0   ;
0   -90   0   40   50   0   0   0   0   0   0   0   0   0   ;
0   0   -40   20   0   20   0   0   0   0   0   0   0   0   ;
40   40   20   -100   0   0   0   0   0   0   0   0   0   0   ;
0   50   0   0   -120   40   0   0   0   0   0   0   0   0   ;
0   0   20   0   40   -70   0   0   0   10   0   0   0   0   ;
10   0   0   0   0   0   -30   0   0   0   0   20   0   0   ;
0   0   0   0   0   0   0   -30   10   0   0   20   0   0   ;
0   0   0   0   0   0   0   10   -60   0   50   0   20   0   ;
0   0   0   0   0   10   0   0   0   -40   10   0   0   0   ;
0   0   0   0   0   0   0   0   50   10   -120   0   0   40   ;
0   0   0   0   0   0   20   20   0   0   0   -50   10   0   ;
0   0   0   0   0   0   0   0   20   0   0   10   -30   0   ;
0   0   0   0   0   0   0   0   0   0   40   0   0   -70   ];

Ywbusdom = [...
-31.5895   4.50587   0   0   10   0   0   0   0   ;
4.50587   -26.1645   0   0   0   0   20   0   0   ;
0   0   -30   10   0   0   20   0   0   ;
0   0   10   -60   0   50   0   20   0   ;
10   0   0   0   -40   10   0   0   0   ;
0   0   0   50   10   -120   0   0   40   ;
0   20   20   0   0   0   -50   10   0   ;
0   0   0   20   0   0   10   -30   0   ;
0   0   0   0   0   40   0   0   -70   ];


Ybusdom = complex(0, Ybusdom);
Ywbusdom = complex(0, Ywbusdom);

Jinjdom = [...
	complex(0,0)	;
	complex(150,-60)	;
	complex(0,0)	;
	complex(140,-60)	;
	complex(-80,-3270)	;
	complex(0,0)	;
	complex(-50,20)	;
	complex(-80,30)	;
	complex(180,-70)	;
	complex(-30,10)	;
	complex(-60,20)	;
	complex(-50,20)	;
	complex(-80,30)	;
	complex(-50,-3280)	];

Jwinjdom = [...
	complex(118.3,-1929.91)	;
	complex(-2.39806,-182.516)	;
	complex(-80,30)	;
	complex(180,-70)	;
	complex(-30,10)	;
	complex(-60,20)	;
	complex(-50,20)	;
	complex(-80,30)	;
	complex(-50,-3280)	];

Yslack = [...
0  0  0  0  30  0  0  0  0  0  0  0  0  30  ];

Ywslack = [...
17.0836  1.6586  0  0  0  0  0  0  30  ];
Yslack = complex(0, Yslack);
Ywslack = complex(0, Ywslack);

Vorig = linsolve(Ybusdom, Jinjdom);
Vw = linsolve(Ywbusdom, Jwinjdom);

slack = 1;
islack_w = Ywslack*Vw + Ywbusdom(slack,slack) * 110;
slack = 6;
islack_orig = Yslack*Vorig + Ybusdom(slack,slack) * 110;
length = 14;

Verr = (abs(Vorig(slack:length))-abs(Vw))./abs(Vorig(slack:length));
Fierr = (angle(Vorig(slack:length))-angle(Vw));
islack_err = (abs(islack_orig) - abs(islack_w))./abs(islack_orig);
islack_ang_err = (angle(islack_orig) - angle(islack_w));
