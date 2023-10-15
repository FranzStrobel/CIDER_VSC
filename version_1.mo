package VSC_FZ

  package Testes
    model Clark3p
      //amplitude invariante
      Modelica.Blocks.Interfaces.RealInput V_abc_D[3] annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-146, -50}, {-106, -10}}, rotation = 0), iconTransformation(origin = {6, 30}, extent = {{-146, -50}, {-106, -10}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealOutput V_ab_D[2] annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{100, -34}, {120, -14}}, rotation = 0), iconTransformation(origin = {0, 24}, extent = {{100, -34}, {120, -14}}, rotation = 0)));
      //
    equation
      V_ab_D[1] = (2/3)*V_abc_D[1] - (1/3)*V_abc_D[2] - (1/3)*V_abc_D[3];
      V_ab_D[2] = 1/sqrt(3)*(V_abc_D[2] - V_abc_D[3]);
      annotation(
        Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}})));
    end Clark3p;

    model VSC
      Modelica.Electrical.Analog.Sources.SignalVoltage Va annotation(
        Placement(visible = true, transformation(origin = {58, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
      Modelica.Electrical.Analog.Sources.SignalVoltage Vb annotation(
        Placement(visible = true, transformation(origin = {60, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
      Modelica.Electrical.Analog.Sources.SignalVoltage Vc annotation(
        Placement(visible = true, transformation(origin = {58, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
      Modelica.Blocks.Math.Product product annotation(
        Placement(visible = true, transformation(origin = {-46, 54}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Product product1 annotation(
        Placement(visible = true, transformation(origin = {-46, -4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Product product2 annotation(
        Placement(visible = true, transformation(origin = {-46, -66}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Interfaces.Pin Vta annotation(
        Placement(visible = true, transformation(origin = {124, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {110, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Interfaces.Pin Vtb annotation(
        Placement(visible = true, transformation(origin = {124, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {110, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Interfaces.Pin Vtc annotation(
        Placement(visible = true, transformation(origin = {124, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {110, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealInput ma annotation(
        Placement(visible = true, transformation(origin = {-133, 59}, extent = {{-13, -13}, {13, 13}}, rotation = 0), iconTransformation(origin = {-60, -114}, extent = {{-14, -14}, {14, 14}}, rotation = 90)));
      Modelica.Blocks.Interfaces.RealInput mb annotation(
        Placement(visible = true, transformation(origin = {-133, 1}, extent = {{-13, -13}, {13, 13}}, rotation = 0), iconTransformation(origin = {0, -114}, extent = {{-14, -14}, {14, 14}}, rotation = 90)));
      Modelica.Blocks.Interfaces.RealInput mc annotation(
        Placement(visible = true, transformation(origin = {-133, -59}, extent = {{-13, -13}, {13, 13}}, rotation = 0), iconTransformation(origin = {60, -114}, extent = {{-14, -14}, {14, 14}}, rotation = 90)));
      Modelica.Electrical.Analog.Interfaces.PositivePin pin_p annotation(
        Placement(visible = true, transformation(origin = {-120, 164}, extent = {{-10, -10}, {10, 10}}, rotation = 180), iconTransformation(origin = {-110, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Interfaces.NegativePin pin_n annotation(
        Placement(visible = true, transformation(origin = {-18, 164}, extent = {{-10, -10}, {10, 10}}, rotation = 180), iconTransformation(origin = {-110, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.VoltageSensor voltageSensor annotation(
        Placement(visible = true, transformation(origin = {-74, 164}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Sources.SignalCurrent I annotation(
        Placement(visible = true, transformation(origin = {-54, 202}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
      Modelica.Electrical.Analog.Sensors.CurrentSensor currentSensor annotation(
        Placement(visible = true, transformation(origin = {92, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.CurrentSensor currentSensor1 annotation(
        Placement(visible = true, transformation(origin = {94, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.CurrentSensor currentSensor2 annotation(
        Placement(visible = true, transformation(origin = {94, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Product product3 annotation(
        Placement(visible = true, transformation(origin = {10, 138}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Product product4 annotation(
        Placement(visible = true, transformation(origin = {24, 110}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Product product5 annotation(
        Placement(visible = true, transformation(origin = {32, 82}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression(y = currentSensor.i) annotation(
        Placement(visible = true, transformation(origin = {-37, 133}, extent = {{-19, -11}, {19, 11}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression1(y = currentSensor1.i) annotation(
        Placement(visible = true, transformation(origin = {-25, 104}, extent = {{-19, -10}, {19, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression2(y = currentSensor2.i) annotation(
        Placement(visible = true, transformation(origin = {-13, 76}, extent = {{-19, -10}, {19, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Add3 add3 annotation(
        Placement(visible = true, transformation(origin = {70, 128}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Division division annotation(
        Placement(visible = true, transformation(origin = {106, 122}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant const(k = 2) annotation(
        Placement(visible = true, transformation(origin = {72, 98}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.CurrentSensor currentSensor3 annotation(
        Placement(visible = true, transformation(origin = {-94, 202}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
      Modelica.Blocks.Math.Division division1 annotation(
        Placement(visible = true, transformation(origin = {-104, 104}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant1(k = 2) annotation(
        Placement(visible = true, transformation(origin = {-144, 98}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(Va.v, product.y) annotation(
        Line(points = {{58, 48}, {58, 38}, {-22, 38}, {-22, 54}, {-34, 54}}, color = {0, 0, 127}));
      connect(Vb.v, product1.y) annotation(
        Line(points = {{60, -12}, {60, -24}, {-22, -24}, {-22, -4}, {-34, -4}}, color = {0, 0, 127}));
      connect(Vc.v, product2.y) annotation(
        Line(points = {{58, -72}, {58, -76}, {-22, -76}, {-22, -66}, {-34, -66}}, color = {0, 0, 127}));
      connect(ma, product.u1) annotation(
        Line(points = {{-133, 59}, {-58, 59}, {-58, 60}}, color = {0, 0, 127}));
      connect(mb, product1.u1) annotation(
        Line(points = {{-133, 1}, {-58, 1}, {-58, 2}}, color = {0, 0, 127}));
      connect(mc, product2.u1) annotation(
        Line(points = {{-133, -59}, {-58, -59}, {-58, -60}}, color = {0, 0, 127}));
      connect(voltageSensor.p, pin_p) annotation(
        Line(points = {{-64, 164}, {-92, 164}, {-92, 164}, {-120, 164}}, color = {0, 0, 255}));
      connect(voltageSensor.n, pin_n) annotation(
        Line(points = {{-84, 164}, {-18, 164}}, color = {0, 0, 255}));
      connect(Va.p, currentSensor.p) annotation(
        Line(points = {{68, 60}, {82, 60}}, color = {0, 0, 255}));
      connect(currentSensor.n, Vta) annotation(
        Line(points = {{102, 60}, {124, 60}}, color = {0, 0, 255}));
      connect(Vb.p, currentSensor1.p) annotation(
        Line(points = {{70, 0}, {84, 0}}, color = {0, 0, 255}));
      connect(Vc.p, currentSensor2.p) annotation(
        Line(points = {{68, -60}, {84, -60}}, color = {0, 0, 255}));
      connect(currentSensor2.n, Vtc) annotation(
        Line(points = {{104, -60}, {124, -60}}, color = {0, 0, 255}));
      connect(currentSensor1.n, Vtb) annotation(
        Line(points = {{104, 0}, {124, 0}}, color = {0, 0, 255}));
      connect(I.p, pin_n) annotation(
        Line(points = {{-44, 202}, {-18, 202}, {-18, 164}}, color = {0, 0, 255}));
      connect(product3.u1, product.u1) annotation(
        Line(points = {{-2, 144}, {-58, 144}, {-58, 60}}, color = {0, 0, 127}));
      connect(realExpression.y, product3.u2) annotation(
        Line(points = {{-16, 133}, {-9, 133}, {-9, 132}, {-2, 132}}, color = {0, 0, 127}));
      connect(realExpression1.y, product4.u2) annotation(
        Line(points = {{-4, 104}, {12, 104}}, color = {0, 0, 127}));
      connect(realExpression2.y, product5.u2) annotation(
        Line(points = {{8, 76}, {20, 76}}, color = {0, 0, 127}));
      connect(product4.u1, product1.u1) annotation(
        Line(points = {{12, 116}, {-66, 116}, {-66, 2}, {-58, 2}}, color = {0, 0, 127}));
      connect(product5.u1, product2.u1) annotation(
        Line(points = {{20, 88}, {-66, 88}, {-66, -60}, {-58, -60}}, color = {0, 0, 127}));
      connect(product3.y, add3.u1) annotation(
        Line(points = {{22, 138}, {58, 138}, {58, 136}}, color = {0, 0, 127}));
      connect(product4.y, add3.u2) annotation(
        Line(points = {{36, 110}, {44, 110}, {44, 128}, {58, 128}}, color = {0, 0, 127}));
      connect(product5.y, add3.u3) annotation(
        Line(points = {{44, 82}, {52, 82}, {52, 120}, {58, 120}}, color = {0, 0, 127}));
      connect(const.y, division.u2) annotation(
        Line(points = {{84, 98}, {86, 98}, {86, 116}, {94, 116}}, color = {0, 0, 127}));
      connect(add3.y, division.u1) annotation(
        Line(points = {{82, 128}, {94, 128}}, color = {0, 0, 127}));
      connect(division.y, I.i) annotation(
        Line(points = {{118, 122}, {136, 122}, {136, 182}, {-54, 182}, {-54, 190}}, color = {0, 0, 127}));
      connect(I.n, currentSensor3.p) annotation(
        Line(points = {{-64, 202}, {-84, 202}}, color = {0, 0, 255}));
      connect(currentSensor3.n, pin_p) annotation(
        Line(points = {{-104, 202}, {-138, 202}, {-138, 164}, {-120, 164}}, color = {0, 0, 255}));
      connect(constant1.y, division1.u2) annotation(
        Line(points = {{-133, 98}, {-116, 98}}, color = {0, 0, 127}));
      connect(division1.y, product.u2) annotation(
        Line(points = {{-93, 104}, {-90, 104}, {-90, 48}, {-58, 48}}, color = {0, 0, 127}));
      connect(division1.y, product1.u2) annotation(
        Line(points = {{-93, 104}, {-90, 104}, {-90, -10}, {-58, -10}}, color = {0, 0, 127}));
      connect(division1.y, product2.u2) annotation(
        Line(points = {{-93, 104}, {-90, 104}, {-90, -72}, {-58, -72}}, color = {0, 0, 127}));
      connect(voltageSensor.v, division1.u1) annotation(
        Line(points = {{-74, 154}, {-74, 124}, {-128, 124}, {-128, 110}, {-116, 110}}, color = {0, 0, 127}));
      connect(Va.n, Vb.n) annotation(
        Line(points = {{48, 60}, {20, 60}, {20, 0}, {50, 0}}, color = {0, 0, 255}));
      connect(Vb.n, Vc.n) annotation(
        Line(points = {{50, 0}, {20, 0}, {20, -60}, {48, -60}}, color = {0, 0, 255}));
      annotation(
        Diagram(coordinateSystem(extent = {{-160, 220}, {140, -80}})));
    end VSC;

    model OuterPQ
      Modelica.Blocks.Interfaces.RealInput vd annotation(
        Placement(visible = true, transformation(origin = {46, -4}, extent = {{-186, 54}, {-146, 94}}, rotation = 0), iconTransformation(origin = {-112, 60}, extent = {{-12, -12}, {12, 12}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealInput vq annotation(
        Placement(visible = true, transformation(origin = {46, -48}, extent = {{-186, 54}, {-146, 94}}, rotation = 0), iconTransformation(origin = {-112, 20}, extent = {{-12, -12}, {12, 12}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealInput pref annotation(
        Placement(visible = true, transformation(origin = {46, -102}, extent = {{-186, 54}, {-146, 94}}, rotation = 0), iconTransformation(origin = {-112, -20}, extent = {{-12, -12}, {12, 12}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealInput qref annotation(
        Placement(visible = true, transformation(origin = {46, -144}, extent = {{-186, 54}, {-146, 94}}, rotation = 0), iconTransformation(origin = {-112, -60}, extent = {{-12, -12}, {12, 12}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealOutput idref annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{100, 46}, {120, 66}}, rotation = 0), iconTransformation(origin = {0, -6}, extent = {{100, 46}, {120, 66}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealOutput iqref annotation(
        Placement(visible = true, transformation(origin = {0, -76}, extent = {{100, 46}, {120, 66}}, rotation = 0), iconTransformation(origin = {0, -106}, extent = {{100, 46}, {120, 66}}, rotation = 0)));
    equation
      idref = (2/3)*(vd*pref + vq*qref)/(vd^2 + vq^2);
      iqref = -(2/3)*(vd*qref - vq*pref)/(vd^2 + vq^2);
    end OuterPQ;

    model InnerControl1
      parameter Real Ictrl_KP;
      parameter Real Ictrl_KI;
      parameter Real Lac_eq_pu;
      parameter Real initvd;
      parameter Real initvq;
      Modelica.Blocks.Interfaces.RealInput Id_ref annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-142, 40}, {-102, 80}}, rotation = 0), iconTransformation(origin = {-110, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealInput Id annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-144, -2}, {-104, 38}}, rotation = 0), iconTransformation(origin = {-110, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealInput Iqref annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-144, -80}, {-104, -40}}, rotation = 0), iconTransformation(origin = {-110, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant const(k = Parameters.lw) annotation(
        Placement(transformation(extent = {{-88, -8}, {-74, 6}})));
      Modelica.Blocks.Math.Product product annotation(
        Placement(transformation(extent = {{-46, 14}, {-34, 26}})));
      Modelica.Blocks.Math.Product product1 annotation(
        Placement(transformation(extent = {{-40, -24}, {-28, -12}})));
      Modelica.Blocks.Math.Add add(k2 = -1) annotation(
        Placement(transformation(extent = {{-74, 50}, {-54, 70}})));
      Modelica.Blocks.Math.Add add1(k1 = -1, k2 = +1) annotation(
        Placement(transformation(extent = {{-78, -56}, {-58, -36}})));
      Modelica.Blocks.Interfaces.RealInput Iq annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-144, -42}, {-104, -2}}, rotation = 0), iconTransformation(origin = {-110, -26}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Add3 add3_1(k2 = +1, k3 = -1) annotation(
        Placement(transformation(extent = {{18, 46}, {38, 66}})));
      Modelica.Blocks.Interfaces.RealInput vd annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-142, 70}, {-102, 110}}, rotation = 0), iconTransformation(origin = {-110, 90}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealInput vq annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-146, -114}, {-106, -74}}, rotation = 0), iconTransformation(origin = {-110, -90}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Add3 add3_2(k1 = +1, k2 = +1) annotation(
        Placement(transformation(extent = {{20, -56}, {40, -36}})));
      Modelica.Blocks.Interfaces.RealOutput vdref annotation(
        Placement(transformation(extent = {{100, 46}, {120, 66}})));
      Modelica.Blocks.Interfaces.RealOutput vqref annotation(
        Placement(transformation(extent = {{100, -56}, {120, -36}})));
      PI1 pi1 annotation(
        Placement(visible = true, transformation(origin = {-22, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      PI1 pi11 annotation(
        Placement(visible = true, transformation(origin = {-26, -46}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Add add2 annotation(
        Placement(visible = true, transformation(origin = {70, 50}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
      Modelica.Blocks.Math.Product product2 annotation(
        Placement(visible = true, transformation(origin = {84, 8}, extent = {{-46, 14}, {-34, 26}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant1(k = Parameters.r) annotation(
        Placement(visible = true, transformation(origin = {100, 8}, extent = {{-88, -8}, {-74, 6}}, rotation = 0)));
      Modelica.Blocks.Math.Product product3 annotation(
        Placement(visible = true, transformation(origin = {82, -36}, extent = {{-46, 14}, {-34, 26}}, rotation = 0)));
      Modelica.Blocks.Math.Add add3 annotation(
        Placement(visible = true, transformation(origin = {74, -48}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
    equation
      connect(const.y, product.u2) annotation(
        Line(points = {{-73.3, -1}, {-56, -1}, {-56, 16.4}, {-47.2, 16.4}}, color = {0, 0, 127}));
      connect(const.y, product1.u1) annotation(
        Line(points = {{-73.3, -1}, {-56, -1}, {-56, -14.4}, {-41.2, -14.4}}, color = {0, 0, 127}));
      connect(Id, add.u2) annotation(
        Line(points = {{-124, 18}, {-90, 18}, {-90, 54}, {-76, 54}}, color = {0, 0, 127}));
      connect(Id_ref, add.u1) annotation(
        Line(points = {{-122, 60}, {-102, 60}, {-76, 60}, {-76, 66}}, color = {0, 0, 127}));
      connect(Id, product.u1) annotation(
        Line(points = {{-124, 18}, {-82, 18}, {-82, 23.6}, {-47.2, 23.6}}, color = {0, 0, 127}));
      connect(Iq, add1.u1) annotation(
        Line(points = {{-124, -22}, {-88, -22}, {-88, -40}, {-80, -40}}, color = {0, 0, 127}));
      connect(Iqref, add1.u2) annotation(
        Line(points = {{-124, -60}, {-106, -60}, {-86, -60}, {-86, -52}, {-80, -52}}, color = {0, 0, 127}));
      connect(Iq, product1.u2) annotation(
        Line(points = {{-124, -22}, {-41.2, -22}, {-41.2, -21.6}}, color = {0, 0, 127}));
      connect(product1.y, add3_1.u3) annotation(
        Line(points = {{-27.4, -18}, {-12, -18}, {6, -18}, {6, 48}, {16, 48}}, color = {0, 0, 127}));
      connect(product.y, add3_2.u1) annotation(
        Line(points = {{-33.4, 20}, {-6, 20}, {-6, -38}, {18, -38}}, color = {0, 0, 127}));
      connect(add.y, pi1.u) annotation(
        Line(points = {{-52, 60}, {-34, 60}}, color = {0, 0, 127}));
      connect(pi1.y, add3_1.u2) annotation(
        Line(points = {{-10, 60}, {6, 60}, {6, 56}, {16, 56}}, color = {0, 0, 127}));
      connect(add1.y, pi11.u) annotation(
        Line(points = {{-56, -46}, {-38, -46}}, color = {0, 0, 127}));
      connect(pi11.y, add3_2.u2) annotation(
        Line(points = {{-14, -46}, {18, -46}}, color = {0, 0, 127}));
      connect(product2.u1, Id) annotation(
        Line(points = {{37, 32}, {-90, 32}, {-90, 18}, {-124, 18}}, color = {0, 0, 127}));
      connect(add3_1.y, add2.u1) annotation(
        Line(points = {{40, 56}, {54, 56}, {54, 54}, {60, 54}}, color = {0, 0, 127}));
      connect(product2.y, add2.u2) annotation(
        Line(points = {{50, 28}, {54, 28}, {54, 46}, {60, 46}}, color = {0, 0, 127}));
      connect(product2.u2, constant1.y) annotation(
        Line(points = {{36, 24}, {32, 24}, {32, 7}, {27, 7}}, color = {0, 0, 127}));
      connect(add2.y, vdref) annotation(
        Line(points = {{78, 50}, {86, 50}, {86, 56}, {110, 56}}, color = {0, 0, 127}));
      connect(product3.u1, constant1.y) annotation(
        Line(points = {{35, -12}, {32, -12}, {32, 8}, {26, 8}}, color = {0, 0, 127}));
      connect(product3.u2, Iq) annotation(
        Line(points = {{34, -20}, {6, -20}, {6, -28}, {-54, -28}, {-54, -22}, {-124, -22}}, color = {0, 0, 127}));
      connect(product3.y, add3.u1) annotation(
        Line(points = {{48, -16}, {56, -16}, {56, -44}, {64, -44}}, color = {0, 0, 127}));
      connect(add3.u2, add3_2.y) annotation(
        Line(points = {{64, -52}, {50, -52}, {50, -46}, {42, -46}}, color = {0, 0, 127}));
      connect(add3.y, vqref) annotation(
        Line(points = {{82, -48}, {93, -48}, {93, -46}, {110, -46}}, color = {0, 0, 127}));
      connect(vd, add3_1.u1) annotation(
        Line(points = {{-122, 90}, {6, 90}, {6, 64}, {16, 64}}, color = {0, 0, 127}));
      connect(vq, add3_2.u3) annotation(
        Line(points = {{-126, -94}, {2, -94}, {2, -54}, {18, -54}}, color = {0, 0, 127}));
      annotation(
        Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}})));
    end InnerControl1;

    model PI1
      parameter Real init_value;
      parameter Real max;
      parameter Real min;
      parameter Real KI;
      parameter Real KP;
      parameter Real time_step = 0.00001;
      Modelica.Blocks.Continuous.LimIntegrator limIntegrator(initType = Modelica.Blocks.Types.Init.InitialOutput, k = 1, limitsAtInit = true, outMax = 1e5, y(start = 0), y_start = init_value) annotation(
        Placement(transformation(extent = {{-14, -10}, {6, 10}})));
      Modelica.Blocks.Interfaces.RealInput u annotation(
        Placement(transformation(extent = {{-142, -20}, {-102, 20}})));
      Modelica.Blocks.Interfaces.RealOutput y annotation(
        Placement(transformation(extent = {{100, -10}, {120, 10}})));
      Modelica.Blocks.Math.Add add(k2 = +1) annotation(
        Placement(transformation(extent = {{-70, -10}, {-50, 10}})));
      Modelica.Blocks.Math.Add add1(k1 = +1, k2 = -1) annotation(
        Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 180, origin = {42, -38})));
      Modelica.Blocks.Math.Add add2(k1 = KP, k2 = KI) annotation(
        Placement(transformation(extent = {{26, 18}, {42, 34}})));
      Modelica.Blocks.Nonlinear.Limiter limiter(limitsAtInit = true, uMax = max, uMin = -max) annotation(
        Placement(transformation(extent = {{58, 20}, {74, 36}})));
    equation
      connect(add.y, limIntegrator.u) annotation(
        Line(points = {{-49, 0}, {-16, 0}}, color = {0, 0, 127}));
      connect(u, add.u1) annotation(
        Line(points = {{-122, 0}, {-106, 0}, {-92, 0}, {-92, 6}, {-72, 6}}, color = {0, 0, 127}));
      connect(add2.u1, add.u1) annotation(
        Line(points = {{24.4, 30.8}, {-92, 30.8}, {-92, 6}, {-72, 6}}, color = {0, 0, 127}));
      connect(limIntegrator.y, add2.u2) annotation(
        Line(points = {{7, 0}, {16, 0}, {16, 21.2}, {24.4, 21.2}}, color = {0, 0, 127}));
      connect(add2.y, limiter.u) annotation(
        Line(points = {{42.8, 26}, {48, 26}, {48, 28}, {56.4, 28}}, color = {0, 0, 127}));
      connect(limiter.y, y) annotation(
        Line(points = {{74.8, 28}, {88, 28}, {88, 0}, {110, 0}}, color = {0, 0, 127}));
      connect(add1.u2, limiter.u) annotation(
        Line(points = {{54, -32}, {58, -32}, {58, -10}, {48, -10}, {48, 28}, {56.4, 28}}, color = {0, 0, 127}));
      connect(add1.u1, y) annotation(
        Line(points = {{54, -44}, {92, -44}, {92, 0}, {110, 0}}, color = {0, 0, 127}));
      connect(add1.y, add.u2) annotation(
        Line(points = {{32, -38}, {-82, -38}, {-82, -6}, {-72, -6}}, color = {0, 0, 127}));
      annotation(
        Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}})));
    end PI1;

    model park_Vinv " park_I block"
      Modelica.Blocks.Interfaces.RealInput Vd annotation(
        Placement(transformation(extent = {{-140, 60}, {-100, 100}})));
      Modelica.Blocks.Interfaces.RealInput Vq annotation(
        Placement(transformation(extent = {{-140, -100}, {-100, -60}})));
      Modelica.Blocks.Interfaces.RealOutput Va annotation(
        Placement(transformation(extent = {{100, 68}, {120, 88}})));
      Modelica.Blocks.Interfaces.RealOutput Vb annotation(
        Placement(transformation(extent = {{100, -88}, {120, -68}})));
      Modelica.Blocks.Interfaces.RealInput theta annotation(
        Placement(transformation(extent = {{20, -20}, {-20, 20}}, rotation = -90, origin = {0, -120})));
    equation
      Va = Vd*sin(theta) + Vq*cos(theta);
      Vb = Vq*sin(theta) - Vd*cos(theta);
      annotation(
        Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics));
    end park_Vinv;

    model InnerControl2
      parameter Real Ictrl_KP;
      parameter Real Ictrl_KI;
      parameter Real Lac_eq_pu;
      parameter Real initvd;
      parameter Real initvq;
      parameter Real r;
      parameter Real lw;
      Modelica.Blocks.Interfaces.RealInput Id_ref annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-142, 40}, {-102, 80}}, rotation = 0), iconTransformation(origin = {-110, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealInput Id annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-144, -2}, {-104, 38}}, rotation = 0), iconTransformation(origin = {-110, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealInput Iqref annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-144, -80}, {-104, -40}}, rotation = 0), iconTransformation(origin = {-110, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant const(k = lw) annotation(
        Placement(transformation(extent = {{-88, -8}, {-74, 6}})));
      Modelica.Blocks.Math.Product product annotation(
        Placement(transformation(extent = {{-46, 14}, {-34, 26}})));
      Modelica.Blocks.Math.Product product1 annotation(
        Placement(transformation(extent = {{-40, -24}, {-28, -12}})));
      Modelica.Blocks.Math.Add add(k2 = -1) annotation(
        Placement(transformation(extent = {{-74, 50}, {-54, 70}})));
      Modelica.Blocks.Math.Add add1(k1 = -1, k2 = +1) annotation(
        Placement(transformation(extent = {{-78, -56}, {-58, -36}})));
      Modelica.Blocks.Interfaces.RealInput Iq annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-144, -42}, {-104, -2}}, rotation = 0), iconTransformation(origin = {-110, -26}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Add3 add3_1(k2 = +1, k3 = -1) annotation(
        Placement(transformation(extent = {{18, 46}, {38, 66}})));
      Modelica.Blocks.Interfaces.RealInput vd annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-142, 70}, {-102, 110}}, rotation = 0), iconTransformation(origin = {-110, 90}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealInput vq annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-146, -114}, {-106, -74}}, rotation = 0), iconTransformation(origin = {-110, -90}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Add3 add3_2(k1 = +1, k2 = +1) annotation(
        Placement(transformation(extent = {{20, -56}, {40, -36}})));
      Modelica.Blocks.Interfaces.RealOutput vdref annotation(
        Placement(transformation(extent = {{100, 46}, {120, 66}})));
      Modelica.Blocks.Interfaces.RealOutput vqref annotation(
        Placement(transformation(extent = {{100, -56}, {120, -36}})));
      PI1 pi1(KI = Ictrl_KI, KP = Ictrl_KP, init_value = 0, max = 1000, min = -1000, time_step = 0.00001) annotation(
        Placement(visible = true, transformation(origin = {-22, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      PI1 pi11(KI = Ictrl_KI, KP = Ictrl_KP, init_value = 0, max = 1000, min = -1000, time_step = 0.00001) annotation(
        Placement(visible = true, transformation(origin = {-26, -46}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Add add2 annotation(
        Placement(visible = true, transformation(origin = {70, 50}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
      Modelica.Blocks.Math.Product product2 annotation(
        Placement(visible = true, transformation(origin = {84, 8}, extent = {{-46, 14}, {-34, 26}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant1(k = r) annotation(
        Placement(visible = true, transformation(origin = {100, 10}, extent = {{-88, -8}, {-74, 6}}, rotation = 0)));
      Modelica.Blocks.Math.Product product3 annotation(
        Placement(visible = true, transformation(origin = {82, -36}, extent = {{-46, 14}, {-34, 26}}, rotation = 0)));
      Modelica.Blocks.Math.Add add3 annotation(
        Placement(visible = true, transformation(origin = {74, -48}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
    equation
      connect(const.y, product.u2) annotation(
        Line(points = {{-73.3, -1}, {-56, -1}, {-56, 16.4}, {-47.2, 16.4}}, color = {0, 0, 127}));
      connect(const.y, product1.u1) annotation(
        Line(points = {{-73.3, -1}, {-56, -1}, {-56, -14.4}, {-41.2, -14.4}}, color = {0, 0, 127}));
      connect(Id, add.u2) annotation(
        Line(points = {{-124, 18}, {-90, 18}, {-90, 54}, {-76, 54}}, color = {0, 0, 127}));
      connect(Id_ref, add.u1) annotation(
        Line(points = {{-122, 60}, {-102, 60}, {-76, 60}, {-76, 66}}, color = {0, 0, 127}));
      connect(Id, product.u1) annotation(
        Line(points = {{-124, 18}, {-82, 18}, {-82, 23.6}, {-47.2, 23.6}}, color = {0, 0, 127}));
      connect(Iq, add1.u1) annotation(
        Line(points = {{-124, -22}, {-88, -22}, {-88, -40}, {-80, -40}}, color = {0, 0, 127}));
      connect(Iqref, add1.u2) annotation(
        Line(points = {{-124, -60}, {-106, -60}, {-86, -60}, {-86, -52}, {-80, -52}}, color = {0, 0, 127}));
      connect(Iq, product1.u2) annotation(
        Line(points = {{-124, -22}, {-41.2, -22}, {-41.2, -21.6}}, color = {0, 0, 127}));
      connect(product1.y, add3_1.u3) annotation(
        Line(points = {{-27.4, -18}, {-12, -18}, {6, -18}, {6, 48}, {16, 48}}, color = {0, 0, 127}));
      connect(product.y, add3_2.u1) annotation(
        Line(points = {{-33.4, 20}, {-6, 20}, {-6, -38}, {18, -38}}, color = {0, 0, 127}));
      connect(add.y, pi1.u) annotation(
        Line(points = {{-52, 60}, {-34, 60}}, color = {0, 0, 127}));
      connect(pi1.y, add3_1.u2) annotation(
        Line(points = {{-10, 60}, {6, 60}, {6, 56}, {16, 56}}, color = {0, 0, 127}));
      connect(add1.y, pi11.u) annotation(
        Line(points = {{-56, -46}, {-38, -46}}, color = {0, 0, 127}));
      connect(pi11.y, add3_2.u2) annotation(
        Line(points = {{-14, -46}, {18, -46}}, color = {0, 0, 127}));
      connect(product2.u1, Id) annotation(
        Line(points = {{37, 32}, {-90, 32}, {-90, 18}, {-124, 18}}, color = {0, 0, 127}));
      connect(add3_1.y, add2.u1) annotation(
        Line(points = {{40, 56}, {54, 56}, {54, 54}, {60, 54}}, color = {0, 0, 127}));
      connect(product2.y, add2.u2) annotation(
        Line(points = {{50, 28}, {54, 28}, {54, 46}, {60, 46}}, color = {0, 0, 127}));
      connect(product2.u2, constant1.y) annotation(
        Line(points = {{36, 24}, {32, 24}, {32, 9}, {27, 9}}, color = {0, 0, 127}));
      connect(add2.y, vdref) annotation(
        Line(points = {{78, 50}, {86, 50}, {86, 56}, {110, 56}}, color = {0, 0, 127}));
      connect(product3.u1, constant1.y) annotation(
        Line(points = {{35, -12}, {32, -12}, {32, 9}, {27, 9}}, color = {0, 0, 127}));
      connect(product3.u2, Iq) annotation(
        Line(points = {{34, -20}, {6, -20}, {6, -28}, {-54, -28}, {-54, -22}, {-124, -22}}, color = {0, 0, 127}));
      connect(product3.y, add3.u1) annotation(
        Line(points = {{48, -16}, {56, -16}, {56, -44}, {64, -44}}, color = {0, 0, 127}));
      connect(add3.u2, add3_2.y) annotation(
        Line(points = {{64, -52}, {50, -52}, {50, -46}, {42, -46}}, color = {0, 0, 127}));
      connect(add3.y, vqref) annotation(
        Line(points = {{82, -48}, {93, -48}, {93, -46}, {110, -46}}, color = {0, 0, 127}));
      connect(vd, add3_1.u1) annotation(
        Line(points = {{-122, 90}, {6, 90}, {6, 64}, {16, 64}}, color = {0, 0, 127}));
      connect(vq, add3_2.u3) annotation(
        Line(points = {{-126, -94}, {2, -94}, {2, -54}, {18, -54}}, color = {0, 0, 127}));
      annotation(
        Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}})));
    end InnerControl2;

    block Park_Fz "Park transformation"
      extends Modelica.Blocks.Icons.Block;
      Modelica.Blocks.Interfaces.RealInput alpha annotation(
        Placement(transformation(extent = {{-140, 20}, {-100, 60}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealInput beta annotation(
        Placement(transformation(extent = {{-140, -60}, {-100, -20}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealOutput d annotation(
        Placement(transformation(extent = {{100, 30}, {120, 50}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealOutput q annotation(
        Placement(transformation(extent = {{100, -50}, {120, -30}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealInput theta annotation(
        Placement(transformation(origin = {0, -120}, extent = {{-20, -20}, {20, 20}}, rotation = 90)));
    equation
      d = alpha*cos(theta) + beta*sin(theta);
      q = -alpha*sin(theta) + beta*cos(theta);
      annotation(
        Diagram(graphics),
        Icon(graphics = {Line(points = {{-75, 0}, {-60, 60}, {-45, 0}, {-30, -60}, {-15, 0}}, color = {0, 0, 255}, smooth = Smooth.Bezier), Line(points = {{-96, -234}}, color = {255, 0, 0}, smooth = Smooth.Bezier), Line(points = {{20, 20}, {80, 20}}, color = {128, 0, 255}), Line(points = {{20, -20}, {80, -20}}, color = {0, 255, 0}), Line(points = {{-60, 0}, {-45, 60}, {-30, 0}, {-15, -60}, {0, 0}}, color = {255, 0, 0}, smooth = Smooth.Bezier)}),
        Documentation(info = "<html>
            <p>
              Perform Park transformation. This transformation translates from the
              static reference frame (alfa-beta) to the synchronous reference
              frame (d-q).
            </p>
          </html>"));
    end Park_Fz;

    block InversePark_fz "Inverse Park transformation"
      extends Modelica.Blocks.Icons.Block;
      Modelica.Blocks.Interfaces.RealInput d annotation(
        Placement(transformation(extent = {{-140, 20}, {-100, 60}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealInput q annotation(
        Placement(transformation(extent = {{-140, -60}, {-100, -20}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealOutput alpha annotation(
        Placement(transformation(extent = {{100, 30}, {120, 50}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealOutput beta annotation(
        Placement(transformation(extent = {{100, -50}, {120, -30}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealInput theta annotation(
        Placement(transformation(origin = {0, -120}, extent = {{-20, -20}, {20, 20}}, rotation = 90)));
    equation
      d = alpha*cos(theta) + beta*sin(theta);
      q = -alpha*sin(theta) + beta*cos(theta);
      annotation(
        Diagram(graphics),
        Icon(graphics = {Line(points = {{0, 0}, {15, 60}, {30, 0}, {45, -60}, {60, 0}}, color = {0, 0, 255}, smooth = Smooth.Bezier), Line(points = {{-96, -234}}, color = {255, 0, 0}, smooth = Smooth.Bezier), Line(points = {{-80, 20}, {-20, 20}}, color = {128, 0, 255}), Line(points = {{-80, -20}, {-20, -20}}, color = {0, 255, 0}), Line(points = {{15, 0}, {30, 60}, {45, 0}, {60, -60}, {75, 0}}, color = {255, 0, 0}, smooth = Smooth.Bezier)}),
        Documentation(info = "<html>
            <p>
              Perform inverse Park transformation. This transformation translates
              from the synchronous reference frame (d-q) to the static reference
              frame (alfa-beta).
            </p>
          </html>"));
    end InversePark_fz;

    block PLL_fz2 "Phase-locked loop"
      extends Modelica.Blocks.Icons.Block;
      parameter Modelica.SIunits.Frequency frequency = 60;
      Modelica.Blocks.Continuous.Integrator integrator annotation(
        Placement(transformation(extent = {{60, -10}, {80, 10}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealInput alpha annotation(
        Placement(visible = true, transformation(origin = {-110, 2}, extent = {{-16, -16}, {16, 16}}, rotation = 0), iconTransformation(origin = {0, 42}, extent = {{-140, -20}, {-100, 20}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealOutput theta annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{100, -10}, {120, 10}}, rotation = 0), iconTransformation(origin = {0, 70}, extent = {{100, -10}, {120, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Add add annotation(
        Placement(transformation(extent = {{30, 48}, {50, 68}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant lineFreq(k = 2*Modelica.Constants.pi*frequency) annotation(
        Placement(visible = true, transformation(origin = {-26, 60}, extent = {{-14, -14}, {14, 14}}, rotation = 0)));
      VSC_FZ.Testes.Park_Fz park_Fz annotation(
        Placement(visible = true, transformation(origin = {-38, 6}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealInput beta annotation(
        Placement(visible = true, transformation(origin = {-109, -25}, extent = {{-17, -17}, {17, 17}}, rotation = 0), iconTransformation(origin = {-2, -46}, extent = {{-140, -20}, {-100, 20}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealOutput d annotation(
        Placement(visible = true, transformation(origin = {0, -52}, extent = {{100, -10}, {120, 10}}, rotation = 0), iconTransformation(origin = {0, 0}, extent = {{100, -10}, {120, 10}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealOutput q annotation(
        Placement(visible = true, transformation(origin = {0, -82}, extent = {{100, -10}, {120, 10}}, rotation = 0), iconTransformation(origin = {0, -64}, extent = {{100, -10}, {120, 10}}, rotation = 0)));
      VSC_FZ.Testes.PI1 pi1(KI = 87.73, KP = 0.9872, init_value = 0, max = 1000, min = -1000) annotation(
        Placement(visible = true, transformation(origin = {6, 24}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(integrator.y, theta) annotation(
        Line(points = {{81, 0}, {110, 0}}, color = {0, 0, 127}));
      connect(add.y, integrator.u) annotation(
        Line(points = {{51, 58}, {54, 58}, {54, 0}, {58, 0}}, color = {0, 0, 127}));
      connect(lineFreq.y, add.u1) annotation(
        Line(points = {{-11, 60}, {9.5, 60}, {9.5, 64}, {28, 64}}, color = {0, 0, 127}));
      connect(park_Fz.alpha, alpha) annotation(
        Line(points = {{-50, 10}, {-81, 10}, {-81, 2}, {-110, 2}}, color = {0, 0, 127}));
      connect(park_Fz.theta, integrator.y) annotation(
        Line(points = {{-38, -6}, {-38, -36}, {90, -36}, {90, 0}, {82, 0}}, color = {0, 0, 127}));
      connect(beta, park_Fz.beta) annotation(
        Line(points = {{-109, -25}, {-70, -25}, {-70, 2}, {-50, 2}}, color = {0, 0, 127}));
      connect(park_Fz.d, d) annotation(
        Line(points = {{-27, 10}, {-12, 10}, {-12, -52}, {110, -52}}, color = {0, 0, 127}));
      connect(park_Fz.q, q) annotation(
        Line(points = {{-27, 2}, {-16, 2}, {-16, -82}, {110, -82}}, color = {0, 0, 127}));
      connect(pi1.y, add.u2) annotation(
        Line(points = {{17, 24}, {17, 52}, {28, 52}}, color = {0, 0, 127}));
      connect(pi1.u, park_Fz.q) annotation(
        Line(points = {{-6, 24}, {-16, 24}, {-16, 2}, {-27, 2}}, color = {0, 0, 127}));
      annotation(
        Diagram(coordinateSystem(extent = {{-140, 80}, {120, -100}})),
        Icon(graphics = {Line(points = {{-70, 0}, {-50, 60}, {-30, 0}, {-10, -60}, {10, 0}, {30, 60}, {50, 0}, {70, -60}, {90, 0}}, color = {0, 0, 255}, smooth = Smooth.Bezier), Line(points = {{-90, 0}, {-64, 60}, {-44, 0}, {-18, -60}, {2, 0}, {22, 60}, {44, 0}, {64, -60}, {88, 0}}, color = {255, 0, 0}, smooth = Smooth.Bezier)}),
        Documentation(info = "<html>
                <p>
                  Phase-locked loop. Given a sinusoidal input, extract the phase.
                </p>
              </html>"));
    end PLL_fz2;

    model PowerControl1
      parameter Real Pctrl_KP;
      parameter Real Pctrl_KI;
      parameter Real Lac_eq_pu;
      parameter Real initvd;
      parameter Real initvq;
      Modelica.Blocks.Interfaces.RealInput P_ref annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-142, 40}, {-102, 80}}, rotation = 0), iconTransformation(origin = {-110, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealInput P annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-144, -2}, {-104, 38}}, rotation = 0), iconTransformation(origin = {-110, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealInput Q_ref annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-144, -80}, {-104, -40}}, rotation = 0), iconTransformation(origin = {-110, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Add add(k2 = -1) annotation(
        Placement(visible = true, transformation(origin = {14, -18}, extent = {{-74, 50}, {-54, 70}}, rotation = 0)));
      Modelica.Blocks.Math.Add add1(k1 = -1, k2 = +1) annotation(
        Placement(visible = true, transformation(origin = {18, 8}, extent = {{-78, -56}, {-58, -36}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealInput Q annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-144, -42}, {-104, -2}}, rotation = 0), iconTransformation(origin = {-110, -26}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealOutput idref annotation(
        Placement(visible = true, transformation(origin = {0, -14}, extent = {{100, 46}, {120, 66}}, rotation = 0), iconTransformation(origin = {0, 0}, extent = {{100, 46}, {120, 66}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealOutput iqref annotation(
        Placement(visible = true, transformation(origin = {0, 8}, extent = {{100, -56}, {120, -36}}, rotation = 0), iconTransformation(origin = {0, 0}, extent = {{100, -56}, {120, -36}}, rotation = 0)));
      VSC_FZ.Testes.PI1 pi1(KI = Pctrl_KI, KP = Pctrl_KP, init_value = 0, max = 1698, min = -1698, time_step = 0.00001) annotation(
        Placement(visible = true, transformation(origin = {24, 42}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      VSC_FZ.Testes.PI1 pi11(KI = -Pctrl_KI, KP = -Pctrl_KP, init_value = 0, max = 1698, min = -1698, time_step = 0.00001) annotation(
        Placement(visible = true, transformation(origin = {24, -38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(P, add.u2) annotation(
        Line(points = {{-124, 18}, {-90, 18}, {-90, 36}, {-62, 36}}, color = {0, 0, 127}));
      connect(P_ref, add.u1) annotation(
        Line(points = {{-122, 60}, {-99, 60}, {-99, 48}, {-62, 48}}, color = {0, 0, 127}));
      connect(Q, add1.u1) annotation(
        Line(points = {{-124, -22}, {-88, -22}, {-88, -32}, {-62, -32}}, color = {0, 0, 127}));
      connect(Q_ref, add1.u2) annotation(
        Line(points = {{-124, -60}, {-86, -60}, {-86, -44}, {-62, -44}}, color = {0, 0, 127}));
      connect(add.y, pi1.u) annotation(
        Line(points = {{-39, 42}, {12, 42}}, color = {0, 0, 127}));
      connect(add1.y, pi11.u) annotation(
        Line(points = {{-39, -38}, {12, -38}}, color = {0, 0, 127}));
      connect(pi1.y, idref) annotation(
        Line(points = {{35, 42}, {110, 42}}, color = {0, 0, 127}));
      connect(pi11.y, iqref) annotation(
        Line(points = {{35, -38}, {110, -38}}, color = {0, 0, 127}));
      annotation(
        Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-140, 80}, {120, -80}})));
    end PowerControl1;

    model PowerCalc
      Modelica.Blocks.Interfaces.RealInput V1 annotation(
        Placement(visible = true, transformation(origin = {2, 28}, extent = {{-142, 40}, {-102, 80}}, rotation = 0), iconTransformation(origin = {-110, 86}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealInput V2 annotation(
        Placement(visible = true, transformation(origin = {2, -2}, extent = {{-142, 40}, {-102, 80}}, rotation = 0), iconTransformation(origin = {-110, 58}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealInput V3 annotation(
        Placement(visible = true, transformation(origin = {2, -34}, extent = {{-142, 40}, {-102, 80}}, rotation = 0), iconTransformation(origin = {-110, 26}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealInput I2 annotation(
        Placement(visible = true, transformation(origin = {2, -108}, extent = {{-142, 40}, {-102, 80}}, rotation = 0), iconTransformation(origin = {-110, -52}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealInput I3 annotation(
        Placement(visible = true, transformation(origin = {2, -142}, extent = {{-142, 40}, {-102, 80}}, rotation = 0), iconTransformation(origin = {-110, -84}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealInput I1 annotation(
        Placement(visible = true, transformation(origin = {2, -74}, extent = {{-142, 40}, {-102, 80}}, rotation = 0), iconTransformation(origin = {-110, -22}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealOutput P annotation(
        Placement(visible = true, transformation(origin = {0, -14}, extent = {{100, 46}, {120, 66}}, rotation = 0), iconTransformation(origin = {0, -10}, extent = {{100, 46}, {120, 66}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealOutput Q annotation(
        Placement(visible = true, transformation(origin = {0, -82}, extent = {{100, 46}, {120, 66}}, rotation = 0), iconTransformation(origin = {0, -86}, extent = {{100, 46}, {120, 66}}, rotation = 0)));
    equation
      P = V1*I1 + V2*I2 + V3*I3;
      Q = ((V1 - V2)*I3 + (V3 - V1)*I2 + (V2 - V3)*I1)/sqrt(3);
    end PowerCalc;

    model y_D_Transformer_ideal
      parameter Real Vac_primary;
      parameter Real Vac_secondary;
      parameter Real trans_ratio = (Vac_primary/Vac_secondary)/sqrt(3);
      Modelica.Electrical.Analog.Ideal.IdealTransformer idealTransformer(n = trans_ratio) annotation(
        Placement(visible = true, transformation(origin = {-2, 52}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Ideal.IdealTransformer idealTransformer1(n = trans_ratio) annotation(
        Placement(visible = true, transformation(origin = {-2, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Ideal.IdealTransformer idealTransformer2(n = trans_ratio) annotation(
        Placement(visible = true, transformation(origin = {0, -42}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Interfaces.NegativePin pin_n annotation(
        Placement(visible = true, transformation(origin = {108, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {110, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Interfaces.NegativePin pin_n1 annotation(
        Placement(visible = true, transformation(origin = {104, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {110, 16}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Interfaces.NegativePin pin_n2 annotation(
        Placement(visible = true, transformation(origin = {108, -34}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {110, -28}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Interfaces.PositivePin pin_p annotation(
        Placement(visible = true, transformation(origin = {-94, 58}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-110, 62}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Interfaces.PositivePin pin_p1 annotation(
        Placement(visible = true, transformation(origin = {-98, 16}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-110, 14}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Interfaces.PositivePin pin_p2 annotation(
        Placement(visible = true, transformation(origin = {-98, -34}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-110, -32}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.VoltageSensor voltageSensor annotation(
        Placement(visible = true, transformation(origin = {-54, 36}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      Modelica.Electrical.Analog.Sensors.VoltageSensor voltageSensor1 annotation(
        Placement(visible = true, transformation(origin = {54, 42}, extent = {{10, -10}, {-10, 10}}, rotation = 90)));
      Modelica.Electrical.Analog.Sensors.CurrentSensor currentSensor annotation(
        Placement(visible = true, transformation(origin = {-48, -6}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Ground ground annotation(
        Placement(visible = true, transformation(origin = {-72, -18}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(idealTransformer.n2, idealTransformer1.p2) annotation(
        Line(points = {{8, 42}, {18, 42}, {18, 12}, {8, 12}}, color = {0, 0, 255}));
      connect(idealTransformer1.n2, idealTransformer2.p2) annotation(
        Line(points = {{8, -8}, {20, -8}, {20, -32}, {10, -32}}, color = {0, 0, 255}));
      connect(idealTransformer.p2, idealTransformer2.n2) annotation(
        Line(points = {{8, 62}, {30, 62}, {30, -52}, {10, -52}}, color = {0, 0, 255}));
      connect(pin_p, idealTransformer.p1) annotation(
        Line(points = {{-94, 58}, {-12, 58}, {-12, 62}}, color = {0, 0, 255}));
      connect(pin_p1, idealTransformer1.p1) annotation(
        Line(points = {{-98, 16}, {-12, 16}, {-12, 12}}, color = {0, 0, 255}));
      connect(pin_p2, idealTransformer2.p1) annotation(
        Line(points = {{-98, -34}, {-10, -34}, {-10, -32}}, color = {0, 0, 255}));
      connect(pin_n, idealTransformer.p2) annotation(
        Line(points = {{108, 60}, {8, 60}, {8, 62}}, color = {0, 0, 255}));
      connect(pin_n1, idealTransformer1.p2) annotation(
        Line(points = {{104, 20}, {8, 20}, {8, 12}}, color = {0, 0, 255}));
      connect(pin_n2, idealTransformer2.p2) annotation(
        Line(points = {{108, -34}, {10, -34}, {10, -32}}, color = {0, 0, 255}));
      connect(voltageSensor.p, idealTransformer.p1) annotation(
        Line(points = {{-54, 46}, {-12, 46}, {-12, 62}}, color = {0, 0, 255}));
      connect(voltageSensor.n, idealTransformer1.p1) annotation(
        Line(points = {{-54, 26}, {-25, 26}, {-25, 12}, {-12, 12}}, color = {0, 0, 255}));
      connect(voltageSensor1.p, idealTransformer.p2) annotation(
        Line(points = {{54, 52}, {54, 62}, {8, 62}}, color = {0, 0, 255}));
      connect(voltageSensor1.n, idealTransformer1.p2) annotation(
        Line(points = {{54, 32}, {50, 32}, {50, 12}, {8, 12}}, color = {0, 0, 255}));
      connect(idealTransformer.n1, currentSensor.p) annotation(
        Line(points = {{-12, 42}, {-32, 42}, {-32, -6}, {-38, -6}}, color = {0, 0, 255}));
      connect(idealTransformer1.n1, currentSensor.p) annotation(
        Line(points = {{-12, -8}, {-38, -8}, {-38, -6}}, color = {0, 0, 255}));
      connect(idealTransformer2.n1, currentSensor.p) annotation(
        Line(points = {{-10, -52}, {-38, -52}, {-38, -6}}, color = {0, 0, 255}));
      connect(ground.p, currentSensor.n) annotation(
        Line(points = {{-72, -8}, {-58, -8}, {-58, -6}}, color = {0, 0, 255}));
      annotation(
        Diagram(coordinateSystem(extent = {{-120, 80}, {120, -60}})));
    end y_D_Transformer_ideal;

    model ART22
      parameter Real Rf = 3.53e-4;
      parameter Real Lf = 2.81e-5;
      parameter Real Rr = 1.47;
      parameter Real Lr = 19.82e-3;
      VSC_FZ.Testes.VSC vsc annotation(
        Placement(visible = true, transformation(origin = {-1, 27}, extent = {{-13, -13}, {13, 13}}, rotation = 0)));
      Modelica.Electrical.Analog.Sources.ConstantVoltage constantVoltage(V = 395.18) annotation(
        Placement(visible = true, transformation(origin = {-86, 18}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      Modelica.Electrical.Analog.Basic.Resistor resistor(R = Rf, alpha = 0) annotation(
        Placement(visible = true, transformation(origin = {87, 34}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Resistor resistor1(R = Rf, alpha = 0) annotation(
        Placement(visible = true, transformation(origin = {90, -26}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Resistor resistor2(R = Rf, alpha = 0) annotation(
        Placement(visible = true, transformation(origin = {92, -86}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Inductor inductor(L = Lf) annotation(
        Placement(visible = true, transformation(origin = {115, 34}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Inductor inductor1(L = Lf) annotation(
        Placement(visible = true, transformation(origin = {118, -26}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      VSC_FZ.Testes.Clark3p clark3p annotation(
        Placement(visible = true, transformation(origin = {-87, 129}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.VoltageSensor vsa annotation(
        Placement(visible = true, transformation(origin = {142, 53}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.VoltageSensor vsb annotation(
        Placement(visible = true, transformation(origin = {141, -6}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.VoltageSensor vscv annotation(
        Placement(visible = true, transformation(origin = {138, -64}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.CurrentSensor ia annotation(
        Placement(visible = true, transformation(origin = {54, 34}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.CurrentSensor ib annotation(
        Placement(visible = true, transformation(origin = {56, -26}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.CurrentSensor ic annotation(
        Placement(visible = true, transformation(origin = {58, -86}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
      VSC_FZ.Testes.Clark3p clark3p1 annotation(
        Placement(visible = true, transformation(origin = {205, 144}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression(y = ia.i) annotation(
        Placement(visible = true, transformation(origin = {171, 144}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression1(y = ib.i) annotation(
        Placement(visible = true, transformation(origin = {171, 128}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression2(y = ic.i) annotation(
        Placement(visible = true, transformation(origin = {171, 112}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression3(y = pLL_fz2.theta) annotation(
        Placement(visible = true, transformation(origin = {215.5, 116}, extent = {{-20.5, -9}, {20.5, 9}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression4(y = vsb.v) annotation(
        Placement(visible = true, transformation(origin = {-129, 123}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression5(y = vscv.v) annotation(
        Placement(visible = true, transformation(origin = {-129, 107}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression6(y = vsa.v) annotation(
        Placement(visible = true, transformation(origin = {-129, 139}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Sources.SineVoltage sineVoltage(V = 11.268e3, freqHz = 60) annotation(
        Placement(visible = true, transformation(origin = {362, 39}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Inductor inductor2(L = Lf) annotation(
        Placement(visible = true, transformation(origin = {120, -86}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Sources.SineVoltage sineVoltage1(V = 11.268e3, freqHz = 60, phase = -2.094395102393195) annotation(
        Placement(visible = true, transformation(origin = {364, -21}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Sources.SineVoltage sineVoltage2(V = 11.268e3, freqHz = 60, phase = 2.094395102393195) annotation(
        Placement(visible = true, transformation(origin = {364, -81}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Ground ground1 annotation(
        Placement(visible = true, transformation(origin = {399, -21}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
      Modelica.Blocks.Sources.RealExpression realExpression7(y = pLL_fz2.d) annotation(
        Placement(visible = true, transformation(origin = {344, 173.5}, extent = {{-17, -9.5}, {17, 9.5}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression8(y = pLL_fz2.q) annotation(
        Placement(visible = true, transformation(origin = {311.5, 98.5}, extent = {{-15.5, -9.5}, {15.5, 9.5}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression9(y = Controle_potencia.idref) annotation(
        Placement(visible = true, transformation(origin = {333, 158.5}, extent = {{-28, -9.5}, {28, 9.5}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression10(y = Correntes.d) annotation(
        Placement(visible = true, transformation(origin = {341, 145}, extent = {{-19, -10}, {19, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression11(y = Correntes.q) annotation(
        Placement(visible = true, transformation(origin = {341.5, 126}, extent = {{-19.5, -10}, {19.5, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression12(y = Controle_potencia.iqref) annotation(
        Placement(visible = true, transformation(origin = {335, 113.5}, extent = {{-26, -8.5}, {26, 8.5}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression13(y = pLL_fz2.theta) annotation(
        Placement(visible = true, transformation(origin = {445, 97.5}, extent = {{-23, -11.5}, {23, 11.5}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant2(k = 0) annotation(
        Placement(visible = true, transformation(origin = {489, 98}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression14(y = clarkInv.a) annotation(
        Placement(visible = true, transformation(origin = {-73, -55.5}, extent = {{-17, -8.5}, {17, 8.5}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression15(y = clarkInv.b) annotation(
        Placement(visible = true, transformation(origin = {-75, -83.5}, extent = {{-17, -9.5}, {17, 9.5}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression16(y = clarkInv.c) annotation(
        Placement(visible = true, transformation(origin = {-76, -111.5}, extent = {{-18, -9.5}, {18, 9.5}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.VoltageSensor vdc_2 annotation(
        Placement(visible = true, transformation(origin = {-109, 21}, extent = {{-7, -7}, {7, 7}}, rotation = -90)));
      Modelica.Blocks.Math.Division division annotation(
        Placement(visible = true, transformation(origin = {434, 157}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Math.Division division1 annotation(
        Placement(visible = true, transformation(origin = {432, 121}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Nonlinear.Limiter limiter(limitsAtInit = true, uMax = 1, uMin = -1) annotation(
        Placement(visible = true, transformation(origin = {-28, -58}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
      Modelica.Blocks.Nonlinear.Limiter limiter1(limitsAtInit = true, uMax = 1, uMin = -1) annotation(
        Placement(visible = true, transformation(origin = {-28, -84}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
      Modelica.Blocks.Nonlinear.Limiter limiter2(limitsAtInit = true, uMax = 1, uMin = -1) annotation(
        Placement(visible = true, transformation(origin = {-28, -112}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
      VSC_FZ.Testes.InnerControl2 Controle_corrente(Ictrl_KI = 0.5545, Ictrl_KP = 4.42e-2, lw = 2*3.1416*60*Lf, r = Rf) annotation(
        Placement(visible = true, transformation(origin = {394.5, 144.5}, extent = {{-14.5, -14.5}, {14.5, 14.5}}, rotation = 0)));
      Modelica.Blocks.Sources.Step Pref(height = 915e3, startTime = 0.1) annotation(
        Placement(visible = true, transformation(origin = {38, 139}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Step Qref(height = 0, startTime = 0) annotation(
        Placement(visible = true, transformation(origin = {37, 73}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      VSC_FZ.Testes.InversePark_fz inversePark_fz annotation(
        Placement(visible = true, transformation(origin = {474, 141}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      VSC_FZ.Testes.Park_Fz Correntes annotation(
        Placement(visible = true, transformation(origin = {243, 140}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      VSC_FZ.Testes.PLL_fz2 pLL_fz2(frequency = 60) annotation(
        Placement(visible = true, transformation(origin = {-34, 133}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Sources.ConstantVoltage constantVoltage1(V = 395.18) annotation(
        Placement(visible = true, transformation(origin = {-86, -16}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      Modelica.Electrical.Analog.Basic.Ground ground annotation(
        Placement(visible = true, transformation(origin = {-65, 70}, extent = {{-54, -86}, {-38, -70}}, rotation = 0)));
      Modelica.Blocks.Continuous.Filter filter(analogFilter = Modelica.Blocks.Types.AnalogFilter.Bessel, f_cut = 60, filterType = Modelica.Blocks.Types.FilterType.LowPass) annotation(
        Placement(visible = true, transformation(origin = {406, -18}, extent = {{84, 4}, {104, 24}}, rotation = 0)));
      VSC_FZ.Testes.PowerCalc powerCalc annotation(
        Placement(visible = true, transformation(origin = {465, -11}, extent = {{-14, -14}, {14, 14}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression17(y = vsa.v) annotation(
        Placement(visible = true, transformation(origin = {431, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression24(y = vsb.v) annotation(
        Placement(visible = true, transformation(origin = {431, -3}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression25(y = vscv.v) annotation(
        Placement(visible = true, transformation(origin = {431, -16}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression26(y = ia.i) annotation(
        Placement(visible = true, transformation(origin = {431, -31}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression27(y = ib.i) annotation(
        Placement(visible = true, transformation(origin = {431, -47}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression28(y = ic.i) annotation(
        Placement(visible = true, transformation(origin = {431, -63}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Continuous.Filter filter1(analogFilter = Modelica.Blocks.Types.AnalogFilter.Bessel, f_cut = 60, filterType = Modelica.Blocks.Types.FilterType.LowPass) annotation(
        Placement(visible = true, transformation(origin = {406, -49}, extent = {{84, 4}, {104, 24}}, rotation = 0)));
      VSC_FZ.Testes.PowerControl1 Controle_potencia(Pctrl_KI = 0.168, Pctrl_KP = 5.36e-4) annotation(
        Placement(visible = true, transformation(origin = {108.5, 104.5}, extent = {{-17.5, -17.5}, {17.5, 17.5}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression23(y = powerCalc.P) annotation(
        Placement(visible = true, transformation(origin = {33.5, 111}, extent = {{-20.5, -10}, {20.5, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression29(y = powerCalc.Q) annotation(
        Placement(visible = true, transformation(origin = {33.5, 94}, extent = {{-20.5, -10}, {20.5, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression18(y = vdc_2.v) annotation(
        Placement(visible = true, transformation(origin = {388, 78.5}, extent = {{-23, -11.5}, {23, 11.5}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Resistor resistor3(R = Rr, alpha = 0) annotation(
        Placement(visible = true, transformation(origin = {299, 39}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Resistor resistor4(R = Rr, alpha = 0) annotation(
        Placement(visible = true, transformation(origin = {298, -21}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Resistor resistor5(R = Rr, alpha = 0) annotation(
        Placement(visible = true, transformation(origin = {298, -81}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Ground ground2 annotation(
        Placement(visible = true, transformation(origin = {161.5, -6.5}, extent = {{-6.5, -6.5}, {6.5, 6.5}}, rotation = 90)));
      Modelica.Electrical.Analog.Basic.Ground ground3 annotation(
        Placement(visible = true, transformation(origin = {160.5, 52.5}, extent = {{-6.5, -6.5}, {6.5, 6.5}}, rotation = 90)));
      Modelica.Electrical.Analog.Basic.Ground ground4 annotation(
        Placement(visible = true, transformation(origin = {160.5, -64.5}, extent = {{-6.5, -6.5}, {6.5, 6.5}}, rotation = 90)));
      Modelica.Electrical.Analog.Basic.Inductor inductor3(L = Lr) annotation(
        Placement(visible = true, transformation(origin = {331, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Inductor inductor4(L = Lr) annotation(
        Placement(visible = true, transformation(origin = {330, -21}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Inductor inductor5(L = Lr) annotation(
        Placement(visible = true, transformation(origin = {330, -81}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      VSC_FZ.Testes.y_D_Transformer_ideal y_D_Transformer_ideal(Vac_primary = 440, Vac_secondary = 13800) annotation(
        Placement(visible = true, transformation(origin = {208, -28}, extent = {{-22, -22}, {22, 22}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Ground ground5 annotation(
        Placement(visible = true, transformation(origin = {302.5, 14.5}, extent = {{-6.5, -6.5}, {6.5, 6.5}}, rotation = 90)));
      Modelica.Electrical.Analog.Sensors.VoltageSensor pca annotation(
        Placement(visible = true, transformation(origin = {284, 15}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.VoltageSensor pcc annotation(
        Placement(visible = true, transformation(origin = {282, -106}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Ground ground6 annotation(
        Placement(visible = true, transformation(origin = {300.5, -106.5}, extent = {{-6.5, -6.5}, {6.5, 6.5}}, rotation = 90)));
      Modelica.Electrical.Analog.Basic.Ground ground7 annotation(
        Placement(visible = true, transformation(origin = {305.5, -44.5}, extent = {{-6.5, -6.5}, {6.5, 6.5}}, rotation = 90)));
      Modelica.Electrical.Analog.Sensors.VoltageSensor pcb annotation(
        Placement(visible = true, transformation(origin = {287, -44}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.CurrentSensor currentSensor annotation(
        Placement(visible = true, transformation(origin = {-97, 1}, extent = {{-4, -4}, {4, 4}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Ground ground8 annotation(
        Placement(visible = true, transformation(origin = {92.5, 9.5}, extent = {{-6.5, -6.5}, {6.5, 6.5}}, rotation = 90)));
      Modelica.Electrical.Analog.Sensors.VoltageSensor vta annotation(
        Placement(visible = true, transformation(origin = {74, 10}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.VoltageSensor vtb annotation(
        Placement(visible = true, transformation(origin = {78, -51}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Ground ground9 annotation(
        Placement(visible = true, transformation(origin = {96.5, -51.5}, extent = {{-6.5, -6.5}, {6.5, 6.5}}, rotation = 90)));
      Modelica.Electrical.Analog.Basic.Ground ground10 annotation(
        Placement(visible = true, transformation(origin = {101.5, -115.5}, extent = {{-6.5, -6.5}, {6.5, 6.5}}, rotation = 90)));
      Modelica.Electrical.Analog.Sensors.VoltageSensor vtc annotation(
        Placement(visible = true, transformation(origin = {83, -115}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Math.Division pcd annotation(
        Placement(visible = true, transformation(origin = {213, 210}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression30(y = pcc.v) annotation(
        Placement(visible = true, transformation(origin = {77, 190}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression31(y = pcb.v) annotation(
        Placement(visible = true, transformation(origin = {77, 206}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression32(y = pca.v) annotation(
        Placement(visible = true, transformation(origin = {77, 222}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Division pcq annotation(
        Placement(visible = true, transformation(origin = {212, 186}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant3(k = 11.267e3) annotation(
        Placement(visible = true, transformation(origin = {173, 185}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      VSC_FZ.Testes.Clark3p clark3p3 annotation(
        Placement(visible = true, transformation(origin = {121, 209}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      VSC_FZ.Testes.PLL_fz2 pLL_fz21(frequency = 60) annotation(
        Placement(visible = true, transformation(origin = {174, 213}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Division vsq annotation(
        Placement(visible = true, transformation(origin = {-16, 81}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Math.Division vsd annotation(
        Placement(visible = true, transformation(origin = {-15, 105}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant4(k = 359.258) annotation(
        Placement(visible = true, transformation(origin = {-55, 80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Division Pref_pu annotation(
        Placement(visible = true, transformation(origin = {130, 159}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant5(k = 915e3) annotation(
        Placement(visible = true, transformation(origin = {97, 143}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant6(k = 915e3) annotation(
        Placement(visible = true, transformation(origin = {506, 31}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Division P_pu annotation(
        Placement(visible = true, transformation(origin = {539, 47}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant7(k = 915e3) annotation(
        Placement(visible = true, transformation(origin = {506, -88}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Division Q_pu annotation(
        Placement(visible = true, transformation(origin = {539, -72}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant8(k = 1698) annotation(
        Placement(visible = true, transformation(origin = {263, 173}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Division id_pu annotation(
        Placement(visible = true, transformation(origin = {296, 189}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant9(k = 1698) annotation(
        Placement(visible = true, transformation(origin = {238, 95}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Division iq_pu annotation(
        Placement(visible = true, transformation(origin = {271, 111}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant10(k = 1698) annotation(
        Placement(visible = true, transformation(origin = {356, 205}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Division idref_pu annotation(
        Placement(visible = true, transformation(origin = {389, 221}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression19(y = Controle_potencia.idref) annotation(
        Placement(visible = true, transformation(origin = {326, 232.5}, extent = {{-28, -9.5}, {28, 9.5}}, rotation = 0)));
      Modelica.Blocks.Math.Division iqref_pu annotation(
        Placement(visible = true, transformation(origin = {503, 221}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant11(k = 1698) annotation(
        Placement(visible = true, transformation(origin = {470, 205}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression21(y = Controle_potencia.iqref) annotation(
        Placement(visible = true, transformation(origin = {444, 229.5}, extent = {{-26, -8.5}, {26, 8.5}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant12(k = 915e3) annotation(
        Placement(visible = true, transformation(origin = {116, 259}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Division Qref_pu annotation(
        Placement(visible = true, transformation(origin = {149, 275}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression20(y = Qref.y) annotation(
        Placement(visible = true, transformation(origin = {107, 283.5}, extent = {{-17, -9.5}, {17, 9.5}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.CurrentSensor ipca annotation(
        Placement(visible = true, transformation(origin = {260, 39}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.CurrentSensor ipcb annotation(
        Placement(visible = true, transformation(origin = {262, -21}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.CurrentSensor ipcc annotation(
        Placement(visible = true, transformation(origin = {264, -81}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression39(y = ia.i) annotation(
        Placement(visible = true, transformation(origin = {-118, 245}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression40(y = ib.i) annotation(
        Placement(visible = true, transformation(origin = {-118, 229}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      VSC_FZ.Testes.PLL_fz2 correntes_tot(frequency = 60) annotation(
        Placement(visible = true, transformation(origin = {-18, 243}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant14(k = 1698) annotation(
        Placement(visible = true, transformation(origin = {-39, 190}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Division itot_d annotation(
        Placement(visible = true, transformation(origin = {1, 215}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Math.Division itot_q annotation(
        Placement(visible = true, transformation(origin = {0, 191}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      VSC_FZ.Testes.Clark3p clark3p2 annotation(
        Placement(visible = true, transformation(origin = {-71, 239}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression38(y = ic.i) annotation(
        Placement(visible = true, transformation(origin = {-118, 213}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant1(k = 915e3) annotation(
        Placement(visible = true, transformation(origin = {655, 34}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression22(y = ipcb.i) annotation(
        Placement(visible = true, transformation(origin = {580, -44}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant13(k = 915e3) annotation(
        Placement(visible = true, transformation(origin = {655, -85}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Continuous.Filter filter2(analogFilter = Modelica.Blocks.Types.AnalogFilter.Bessel, f_cut = 60, filterType = Modelica.Blocks.Types.FilterType.LowPass) annotation(
        Placement(visible = true, transformation(origin = {555, -46}, extent = {{84, 4}, {104, 24}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression33(y = ipca.i) annotation(
        Placement(visible = true, transformation(origin = {580, -28}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Division P_pupc annotation(
        Placement(visible = true, transformation(origin = {688, 50}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Continuous.Filter filter3(analogFilter = Modelica.Blocks.Types.AnalogFilter.Bessel, f_cut = 60, filterType = Modelica.Blocks.Types.FilterType.LowPass) annotation(
        Placement(visible = true, transformation(origin = {555, -15}, extent = {{84, 4}, {104, 24}}, rotation = 0)));
      Modelica.Blocks.Math.Division Q_pupc annotation(
        Placement(visible = true, transformation(origin = {688, -69}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression34(y = pcb.v) annotation(
        Placement(visible = true, transformation(origin = {580, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression35(y = ipcc.i) annotation(
        Placement(visible = true, transformation(origin = {580, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression36(y = pcc.v) annotation(
        Placement(visible = true, transformation(origin = {580, -13}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      VSC_FZ.Testes.PowerCalc powerCalc_pcc annotation(
        Placement(visible = true, transformation(origin = {614, -8}, extent = {{-14, -14}, {14, 14}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression37(y = pca.v) annotation(
        Placement(visible = true, transformation(origin = {580, 13}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      VSC_FZ.Testes.PLL_fz2 ipc(frequency = 60) annotation(
        Placement(visible = true, transformation(origin = {611, 276}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression41(y = ipcc.i) annotation(
        Placement(visible = true, transformation(origin = {516, 243}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Division ipd annotation(
        Placement(visible = true, transformation(origin = {677, 243}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression42(y = ipcb.i) annotation(
        Placement(visible = true, transformation(origin = {511, 262}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      VSC_FZ.Testes.Clark3p clark3p4 annotation(
        Placement(visible = true, transformation(origin = {560, 274}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant15(k = 54.13) annotation(
        Placement(visible = true, transformation(origin = {630, 217}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression43(y = ipca.i) annotation(
        Placement(visible = true, transformation(origin = {511, 279}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Division ipq annotation(
        Placement(visible = true, transformation(origin = {669, 218}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
  ClarkInv clarkInv annotation(
        Placement(visible = true, transformation(origin = {533, 139}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(constantVoltage.p, vsc.pin_p) annotation(
        Line(points = {{-86, 28}, {-86, 34}, {-15, 34}, {-15, 33.5}}, color = {0, 0, 255}));
      connect(resistor.n, inductor.p) annotation(
        Line(points = {{97, 34}, {105, 34}}, color = {0, 0, 255}));
      connect(resistor1.n, inductor1.p) annotation(
        Line(points = {{100, -26}, {108, -26}}, color = {0, 0, 255}));
      connect(vsc.Vta, ia.p) annotation(
        Line(points = {{13, 35}, {30, 35}, {30, 34}, {48, 34}}, color = {0, 0, 255}));
      connect(ia.n, resistor.p) annotation(
        Line(points = {{60, 34}, {77, 34}}, color = {0, 0, 255}));
      connect(vsc.Vtb, ib.p) annotation(
        Line(points = {{13, 27}, {40, 27}, {40, -26}, {50, -26}}, color = {0, 0, 255}));
      connect(ib.n, resistor1.p) annotation(
        Line(points = {{62, -26}, {80, -26}}, color = {0, 0, 255}));
      connect(ic.n, resistor2.p) annotation(
        Line(points = {{64, -86}, {82, -86}}, color = {0, 0, 255}));
      connect(ic.p, vsc.Vtc) annotation(
        Line(points = {{52, -86}, {24, -86}, {24, 19}, {13, 19}}, color = {0, 0, 255}));
      connect(realExpression.y, clark3p1.V_abc_D[1]) annotation(
        Line(points = {{182, 144}, {192, 144}}, color = {0, 0, 127}));
      connect(realExpression1.y, clark3p1.V_abc_D[2]) annotation(
        Line(points = {{182, 128}, {184, 128}, {184, 144}, {192, 144}}, color = {0, 0, 127}));
      connect(realExpression2.y, clark3p1.V_abc_D[3]) annotation(
        Line(points = {{182, 112}, {184, 112}, {184, 144}, {192, 144}}, color = {0, 0, 127}));
      connect(realExpression6.y, clark3p.V_abc_D[1]) annotation(
        Line(points = {{-118, 139}, {-104, 139}, {-104, 129}, {-100, 129}}, color = {0, 0, 127}));
      connect(realExpression4.y, clark3p.V_abc_D[2]) annotation(
        Line(points = {{-118, 123}, {-112, 123}, {-112, 129}, {-100, 129}}, color = {0, 0, 127}));
      connect(realExpression5.y, clark3p.V_abc_D[3]) annotation(
        Line(points = {{-118, 107}, {-112, 107}, {-112, 129}, {-100, 129}}, color = {0, 0, 127}));
      connect(inductor2.p, resistor2.n) annotation(
        Line(points = {{110, -86}, {102, -86}}, color = {0, 0, 255}));
      connect(sineVoltage.n, ground1.p) annotation(
        Line(points = {{372, 39}, {389, 39}, {389, -21}}, color = {0, 0, 255}));
      connect(sineVoltage1.n, ground1.p) annotation(
        Line(points = {{374, -21}, {389, -21}}, color = {0, 0, 255}));
      connect(sineVoltage2.n, ground1.p) annotation(
        Line(points = {{374, -81}, {389, -81}, {389, -21}}, color = {0, 0, 255}));
      connect(vsa.p, inductor.n) annotation(
        Line(points = {{135, 53}, {125, 53}, {125, 34}}, color = {0, 0, 255}));
      connect(vsb.p, inductor1.n) annotation(
        Line(points = {{134, -6}, {128, -6}, {128, -26}}, color = {0, 0, 255}));
      connect(vscv.p, inductor2.n) annotation(
        Line(points = {{130, -64}, {130, -86}}, color = {0, 0, 255}));
      connect(vdc_2.p, constantVoltage.p) annotation(
        Line(points = {{-109, 28}, {-87, 28}}, color = {0, 0, 255}));
      connect(vdc_2.n, constantVoltage.n) annotation(
        Line(points = {{-109, 14}, {-109, 8}, {-87, 8}}, color = {0, 0, 255}));
      connect(realExpression14.y, limiter.u) annotation(
        Line(points = {{-54, -55.5}, {-46, -55.5}, {-46, -58}, {-38, -58}}, color = {0, 0, 127}));
      connect(limiter.y, vsc.ma) annotation(
        Line(points = {{-19, -58}, {-8, -58}, {-8, 12}}, color = {0, 0, 127}));
      connect(realExpression15.y, limiter1.u) annotation(
        Line(points = {{-56, -83.5}, {-47, -83.5}, {-47, -84}, {-38, -84}}, color = {0, 0, 127}));
      connect(realExpression16.y, limiter2.u) annotation(
        Line(points = {{-56, -111.5}, {-47, -111.5}, {-47, -112}, {-38, -112}}, color = {0, 0, 127}));
      connect(limiter1.y, vsc.mb) annotation(
        Line(points = {{-20, -84}, {0, -84}, {0, 12}}, color = {0, 0, 127}));
      connect(limiter2.y, vsc.mc) annotation(
        Line(points = {{-20, -112}, {6, -112}, {6, 12}}, color = {0, 0, 127}));
      connect(Controle_corrente.vdref, division.u1) annotation(
        Line(points = {{410, 153}, {410, 161.22}, {424.45, 161.22}}, color = {0, 0, 127}));
      connect(Controle_corrente.vqref, division1.u1) annotation(
        Line(points = {{410, 138}, {410, 125.23}, {422.45, 125.23}}, color = {0, 0, 127}));
      connect(realExpression7.y, Controle_corrente.vd) annotation(
        Line(points = {{362.7, 173.5}, {373.7, 173.5}, {373.7, 158}, {379, 158}}, color = {0, 0, 127}));
      connect(realExpression10.y, Controle_corrente.Id) annotation(
        Line(points = {{361.9, 145}, {370.4, 145}, {370.4, 146}, {379, 146}}, color = {0, 0, 127}));
      connect(realExpression11.y, Controle_corrente.Iq) annotation(
        Line(points = {{362.95, 126}, {369.95, 126}, {369.95, 141}, {379, 141}}, color = {0, 0, 127}));
      connect(realExpression8.y, Controle_corrente.vq) annotation(
        Line(points = {{328.55, 98.5}, {379, 98.5}, {379, 131}}, color = {0, 0, 127}));
      connect(division.y, inversePark_fz.d) annotation(
        Line(points = {{441.7, 157}, {453.7, 157}, {453.7, 145}, {461.7, 145}}, color = {0, 0, 127}));
      connect(division1.y, inversePark_fz.q) annotation(
        Line(points = {{439.7, 121}, {455.7, 121}, {455.7, 137}, {461.7, 137}}, color = {0, 0, 127}));
      connect(inversePark_fz.theta, realExpression13.y) annotation(
        Line(points = {{474, 129}, {474, 97.5}, {470, 97.5}}, color = {0, 0, 127}));
      connect(clark3p1.V_ab_D[1], Correntes.alpha) annotation(
        Line(points = {{216, 144}, {231, 144}}, color = {0, 0, 127}));
      connect(clark3p1.V_ab_D[2], Correntes.beta) annotation(
        Line(points = {{216, 144}, {224, 144}, {224, 136}, {231, 136}}, color = {0, 0, 127}));
      connect(Correntes.theta, realExpression3.y) annotation(
        Line(points = {{243, 128}, {243, 116}, {238, 116}}, color = {0, 0, 127}));
      connect(clark3p.V_ab_D[1], pLL_fz2.alpha) annotation(
        Line(points = {{-76, 129}, {-63.5, 129}, {-63.5, 137}, {-46, 137}}, color = {0, 0, 127}));
      connect(pLL_fz2.beta, clark3p.V_ab_D[2]) annotation(
        Line(points = {{-46.2, 128.4}, {-59.2, 128.4}, {-59.2, 129.4}, {-76.2, 129.4}}, color = {0, 0, 127}));
      connect(constantVoltage1.n, vsc.pin_n) annotation(
        Line(points = {{-86, -26}, {-60, -26}, {-60, 20.5}, {-15, 20.5}}, color = {0, 0, 255}));
      connect(constantVoltage.n, constantVoltage1.p) annotation(
        Line(points = {{-86, 8}, {-86, -6}}, color = {0, 0, 255}));
      connect(realExpression17.y, powerCalc.V1) annotation(
        Line(points = {{442, 10}, {447, 10}, {447, 1}, {450, 1}}, color = {0, 0, 127}));
      connect(realExpression24.y, powerCalc.V2) annotation(
        Line(points = {{442, -3}, {450, -3}}, color = {0, 0, 127}));
      connect(realExpression25.y, powerCalc.V3) annotation(
        Line(points = {{442, -16}, {444, -16}, {444, -7}, {450, -7}}, color = {0, 0, 127}));
      connect(realExpression26.y, powerCalc.I1) annotation(
        Line(points = {{442, -31}, {446, -31}, {446, -14}, {450, -14}}, color = {0, 0, 127}));
      connect(realExpression27.y, powerCalc.I2) annotation(
        Line(points = {{442, -47}, {447, -47}, {447, -18}, {450, -18}}, color = {0, 0, 127}));
      connect(realExpression28.y, powerCalc.I3) annotation(
        Line(points = {{442, -63}, {448, -63}, {448, -23}, {450, -23}}, color = {0, 0, 127}));
      connect(powerCalc.P, filter.u) annotation(
        Line(points = {{480.4, -4.56}, {488.4, -4.56}, {488.4, -3.56}}, color = {0, 0, 127}));
      connect(powerCalc.Q, filter1.u) annotation(
        Line(points = {{480.4, -15.2}, {483.4, -15.2}, {483.4, -35.2}, {488.4, -35.2}}, color = {0, 0, 127}));
      connect(Qref.y, Controle_potencia.Q_ref) annotation(
        Line(points = {{48, 73}, {78, 73}, {78, 94}, {89, 94}}, color = {0, 0, 127}));
      connect(Pref.y, Controle_potencia.P_ref) annotation(
        Line(points = {{49, 139}, {81, 139}, {81, 113}, {89, 113}}, color = {0, 0, 127}));
      connect(realExpression9.y, Controle_corrente.Id_ref) annotation(
        Line(points = {{364, 158.5}, {369.25, 158.5}, {369.25, 152}, {379, 152}}, color = {0, 0, 127}));
      connect(realExpression12.y, Controle_corrente.Iqref) annotation(
        Line(points = {{364, 113.5}, {374.1, 113.5}, {374.1, 136}, {379, 136}}, color = {0, 0, 127}));
      connect(realExpression29.y, Controle_potencia.Q) annotation(
        Line(points = {{56, 94}, {72, 94}, {72, 100}, {89, 100}}, color = {0, 0, 127}));
      connect(realExpression23.y, Controle_potencia.P) annotation(
        Line(points = {{56, 111}, {67, 111}, {67, 106}, {89, 106}}, color = {0, 0, 127}));
      connect(realExpression18.y, division1.u2) annotation(
        Line(points = {{413.3, 78.5}, {413.3, 117}, {424.3, 117}}, color = {0, 0, 127}));
      connect(realExpression18.y, division.u2) annotation(
        Line(points = {{413.3, 78.5}, {413.3, 153}, {426.3, 153}}, color = {0, 0, 127}));
      connect(vsb.n, ground2.p) annotation(
        Line(points = {{148, -6}, {155, -6}}, color = {0, 0, 255}));
      connect(vsa.n, ground3.p) annotation(
        Line(points = {{149, 53}, {154, 53}}, color = {0, 0, 255}));
      connect(vscv.n, ground4.p) annotation(
        Line(points = {{146, -64}, {154, -64}}, color = {0, 0, 255}));
      connect(resistor4.n, inductor4.n) annotation(
        Line(points = {{308, -21}, {320, -21}}, color = {0, 0, 255}));
      connect(resistor5.n, inductor5.n) annotation(
        Line(points = {{308, -81}, {320, -81}}, color = {0, 0, 255}));
      connect(inductor5.p, sineVoltage2.p) annotation(
        Line(points = {{340, -81}, {354, -81}}, color = {0, 0, 255}));
      connect(inductor4.p, sineVoltage1.p) annotation(
        Line(points = {{340, -21}, {354, -21}}, color = {0, 0, 255}));
      connect(resistor3.n, inductor3.p) annotation(
        Line(points = {{309, 39}, {313, 39}, {313, 40}, {321, 40}}, color = {0, 0, 255}));
      connect(inductor3.n, sineVoltage.p) annotation(
        Line(points = {{341, 40}, {352, 40}, {352, 39}}, color = {0, 0, 255}));
      connect(y_D_Transformer_ideal.pin_p, inductor.n) annotation(
        Line(points = {{184, -14}, {180, -14}, {180, 34}, {125, 34}}, color = {0, 0, 255}));
      connect(y_D_Transformer_ideal.pin_p1, inductor1.n) annotation(
        Line(points = {{184, -25}, {128, -25}, {128, -26}}, color = {0, 0, 255}));
      connect(y_D_Transformer_ideal.pin_p2, inductor2.n) annotation(
        Line(points = {{184, -35}, {171, -35}, {171, -86}, {130, -86}}, color = {0, 0, 255}));
      connect(pca.n, ground5.p) annotation(
        Line(points = {{291, 15}, {296, 15}}, color = {0, 0, 255}));
      connect(pcc.n, ground6.p) annotation(
        Line(points = {{289, -106}, {294, -106}}, color = {0, 0, 255}));
      connect(pcb.n, ground7.p) annotation(
        Line(points = {{294, -44}, {299, -44}}, color = {0, 0, 255}));
      connect(pca.p, y_D_Transformer_ideal.pin_n) annotation(
        Line(points = {{277, 15}, {238, 15}, {238, -15}, {232, -15}}, color = {0, 0, 255}));
      connect(pcb.p, y_D_Transformer_ideal.pin_n1) annotation(
        Line(points = {{280, -44}, {248, -44}, {248, -24}, {232, -24}}, color = {0, 0, 255}));
      connect(pcc.p, y_D_Transformer_ideal.pin_n2) annotation(
        Line(points = {{275, -106}, {243, -106}, {243, -34}, {232, -34}}, color = {0, 0, 255}));
      connect(ground.p, currentSensor.p) annotation(
        Line(points = {{-111, 0}, {-101, 0}, {-101, 1}}, color = {0, 0, 255}));
      connect(currentSensor.n, constantVoltage.n) annotation(
        Line(points = {{-93, 1}, {-86, 1}, {-86, 8}}, color = {0, 0, 255}));
      connect(vta.n, ground8.p) annotation(
        Line(points = {{81, 10}, {86, 10}}, color = {0, 0, 255}));
      connect(vtb.n, ground9.p) annotation(
        Line(points = {{85, -51}, {90, -51}}, color = {0, 0, 255}));
      connect(vtc.n, ground10.p) annotation(
        Line(points = {{90, -115}, {95, -115}}, color = {0, 0, 255}));
      connect(vtc.p, resistor2.p) annotation(
        Line(points = {{76, -115}, {70, -115}, {70, -86}, {82, -86}}, color = {0, 0, 255}));
      connect(vtb.p, resistor1.p) annotation(
        Line(points = {{71, -51}, {65, -51}, {65, -26}, {80, -26}}, color = {0, 0, 255}));
      connect(vta.p, resistor.p) annotation(
        Line(points = {{67, 10}, {64, 10}, {64, 34}, {77, 34}}, color = {0, 0, 255}));
      connect(realExpression31.y, clark3p3.V_abc_D[2]) annotation(
        Line(points = {{88, 206}, {94, 206}, {94, 209}, {109, 209}}, color = {0, 0, 127}));
      connect(pLL_fz21.q, pcq.u1) annotation(
        Line(points = {{185, 206.6}, {193, 206.6}, {193, 189.6}, {204, 189.6}}, color = {0, 0, 127}));
      connect(realExpression32.y, clark3p3.V_abc_D[1]) annotation(
        Line(points = {{88, 222}, {102, 222}, {102, 209}, {109, 209}}, color = {0, 0, 127}));
      connect(realExpression30.y, clark3p3.V_abc_D[3]) annotation(
        Line(points = {{88, 190}, {94, 190}, {94, 209}, {109, 209}}, color = {0, 0, 127}));
      connect(pLL_fz21.d, pcd.u1) annotation(
        Line(points = {{185, 213}, {205, 213}, {205, 214}}, color = {0, 0, 127}));
      connect(clark3p3.V_ab_D[1], pLL_fz21.alpha) annotation(
        Line(points = {{132, 209}, {144.5, 209}, {144.5, 217}, {162, 217}}, color = {0, 0, 127}));
      connect(constant3.y, pcq.u2) annotation(
        Line(points = {{184, 185}, {204, 185}, {204, 182}}, color = {0, 0, 127}));
      connect(pLL_fz21.beta, clark3p3.V_ab_D[2]) annotation(
        Line(points = {{161.8, 208.4}, {148.8, 208.4}, {148.8, 209.4}, {131.8, 209.4}}, color = {0, 0, 127}));
      connect(constant3.y, pcd.u2) annotation(
        Line(points = {{184, 185}, {197, 185}, {197, 206}, {205, 206}}, color = {0, 0, 127}));
      connect(constant4.y, vsq.u2) annotation(
        Line(points = {{-44, 80}, {-24, 80}, {-24, 77}}, color = {0, 0, 127}));
      connect(constant4.y, vsd.u2) annotation(
        Line(points = {{-44, 80}, {-31, 80}, {-31, 101}, {-23, 101}}, color = {0, 0, 127}));
      connect(pLL_fz2.d, vsd.u1) annotation(
        Line(points = {{-23, 133}, {-33, 133}, {-33, 109}, {-23, 109}}, color = {0, 0, 127}));
      connect(pLL_fz2.q, vsq.u1) annotation(
        Line(points = {{-23, 127}, {-37, 127}, {-37, 85}, {-24, 85}}, color = {0, 0, 127}));
      connect(Pref_pu.u1, Pref.y) annotation(
        Line(points = {{122, 163}, {81, 163}, {81, 139}, {49, 139}}, color = {0, 0, 127}));
      connect(constant5.y, Pref_pu.u2) annotation(
        Line(points = {{108, 143}, {113, 143}, {113, 155}, {122, 155}}, color = {0, 0, 127}));
      connect(constant6.y, P_pu.u2) annotation(
        Line(points = {{517, 31}, {522, 31}, {522, 43}, {531, 43}}, color = {0, 0, 127}));
      connect(powerCalc.P, P_pu.u1) annotation(
        Line(points = {{480, -5}, {482, -5}, {482, 51}, {531, 51}}, color = {0, 0, 127}));
      connect(constant7.y, Q_pu.u2) annotation(
        Line(points = {{517, -88}, {522, -88}, {522, -76}, {531, -76}}, color = {0, 0, 127}));
      connect(powerCalc.Q, Q_pu.u1) annotation(
        Line(points = {{480, -15}, {484, -15}, {484, -68}, {531, -68}}, color = {0, 0, 127}));
      connect(constant8.y, id_pu.u2) annotation(
        Line(points = {{274, 173}, {279, 173}, {279, 185}, {288, 185}}, color = {0, 0, 127}));
      connect(id_pu.u1, Correntes.d) annotation(
        Line(points = {{288, 193}, {242, 193}, {242, 152}, {263, 152}, {263, 144}, {254, 144}}, color = {0, 0, 127}));
      connect(constant9.y, iq_pu.u2) annotation(
        Line(points = {{249, 95}, {254, 95}, {254, 107}, {263, 107}}, color = {0, 0, 127}));
      connect(iq_pu.u1, Correntes.q) annotation(
        Line(points = {{263, 115}, {257, 115}, {257, 136}, {254, 136}}, color = {0, 0, 127}));
      connect(constant10.y, idref_pu.u2) annotation(
        Line(points = {{367, 205}, {372, 205}, {372, 217}, {381, 217}}, color = {0, 0, 127}));
      connect(realExpression19.y, idref_pu.u1) annotation(
        Line(points = {{357, 233}, {375, 233}, {375, 225}, {381, 225}}, color = {0, 0, 127}));
      connect(constant11.y, iqref_pu.u2) annotation(
        Line(points = {{481, 205}, {486, 205}, {486, 217}, {495, 217}}, color = {0, 0, 127}));
      connect(realExpression21.y, iqref_pu.u1) annotation(
        Line(points = {{473, 229.5}, {484, 229.5}, {484, 225}, {495, 225}}, color = {0, 0, 127}));
      connect(constant12.y, Qref_pu.u2) annotation(
        Line(points = {{127, 259}, {132, 259}, {132, 271}, {141, 271}}, color = {0, 0, 127}));
      connect(realExpression20.y, Qref_pu.u1) annotation(
        Line(points = {{126, 284}, {133, 284}, {133, 279}, {141, 279}}, color = {0, 0, 127}));
      connect(resistor3.p, ipca.n) annotation(
        Line(points = {{289, 39}, {266, 39}}, color = {0, 0, 255}));
      connect(ipca.p, y_D_Transformer_ideal.pin_n) annotation(
        Line(points = {{254, 39}, {238, 39}, {238, -15}, {232, -15}}, color = {0, 0, 255}));
      connect(resistor4.p, ipcb.n) annotation(
        Line(points = {{288, -21}, {268, -21}}, color = {0, 0, 255}));
      connect(ipcb.p, y_D_Transformer_ideal.pin_n1) annotation(
        Line(points = {{256, -21}, {248, -21}, {248, -24}, {232, -24}}, color = {0, 0, 255}));
      connect(resistor5.p, ipcc.n) annotation(
        Line(points = {{288, -81}, {270, -81}}, color = {0, 0, 255}));
      connect(ipcc.p, y_D_Transformer_ideal.pin_n2) annotation(
        Line(points = {{258, -81}, {243, -81}, {243, -34}, {232, -34}}, color = {0, 0, 255}));
      connect(constant14.y, itot_q.u2) annotation(
        Line(points = {{-28, 190}, {-8, 190}, {-8, 187}}, color = {0, 0, 127}));
      connect(realExpression39.y, clark3p2.V_abc_D[1]) annotation(
        Line(points = {{-107, 245}, {-91, 245}, {-91, 239}, {-83, 239}}, color = {0, 0, 127}));
      connect(clark3p2.V_ab_D[1], correntes_tot.alpha) annotation(
        Line(points = {{-60, 239}, {-47.5, 239}, {-47.5, 247}, {-30, 247}}, color = {0, 0, 127}));
      connect(correntes_tot.q, itot_q.u1) annotation(
        Line(points = {{-7, 236.6}, {-21, 236.6}, {-21, 195}, {-8, 195}}, color = {0, 0, 127}));
      connect(constant14.y, itot_d.u2) annotation(
        Line(points = {{-28, 190}, {-15, 190}, {-15, 211}, {-7, 211}}, color = {0, 0, 127}));
      connect(realExpression40.y, clark3p2.V_abc_D[2]) annotation(
        Line(points = {{-107, 229}, {-91, 229}, {-91, 239}, {-83, 239}}, color = {0, 0, 127}));
      connect(correntes_tot.beta, clark3p2.V_ab_D[2]) annotation(
        Line(points = {{-30.2, 238.4}, {-43.2, 238.4}, {-43.2, 239.4}, {-60.2, 239.4}}, color = {0, 0, 127}));
      connect(realExpression38.y, clark3p2.V_abc_D[3]) annotation(
        Line(points = {{-107, 213}, {-91, 213}, {-91, 239}, {-83, 239}}, color = {0, 0, 127}));
      connect(correntes_tot.d, itot_d.u1) annotation(
        Line(points = {{-7, 243}, {-17, 243}, {-17, 219}, {-7, 219}}, color = {0, 0, 127}));
      connect(constant1.y, P_pupc.u2) annotation(
        Line(points = {{666, 34}, {671, 34}, {671, 46}, {680, 46}}, color = {0, 0, 127}));
      connect(realExpression33.y, powerCalc_pcc.I1) annotation(
        Line(points = {{591, -28}, {595, -28}, {595, -11}, {599, -11}}, color = {0, 0, 127}));
      connect(realExpression35.y, powerCalc_pcc.I3) annotation(
        Line(points = {{591, -60}, {597, -60}, {597, -20}, {599, -20}}, color = {0, 0, 127}));
      connect(constant13.y, Q_pupc.u2) annotation(
        Line(points = {{666, -85}, {671, -85}, {671, -73}, {680, -73}}, color = {0, 0, 127}));
      connect(realExpression37.y, powerCalc_pcc.V1) annotation(
        Line(points = {{591, 13}, {596, 13}, {596, 4}, {599, 4}}, color = {0, 0, 127}));
      connect(powerCalc_pcc.Q, Q_pupc.u1) annotation(
        Line(points = {{629.4, -12.2}, {633.4, -12.2}, {633.4, -65.2}, {680.4, -65.2}}, color = {0, 0, 127}));
      connect(realExpression22.y, powerCalc_pcc.I2) annotation(
        Line(points = {{591, -44}, {596, -44}, {596, -15}, {599, -15}}, color = {0, 0, 127}));
      connect(realExpression34.y, powerCalc_pcc.V2) annotation(
        Line(points = {{591, 0}, {599, 0}}, color = {0, 0, 127}));
      connect(realExpression36.y, powerCalc_pcc.V3) annotation(
        Line(points = {{591, -13}, {593, -13}, {593, -4}, {599, -4}}, color = {0, 0, 127}));
      connect(powerCalc_pcc.Q, filter2.u) annotation(
        Line(points = {{629.4, -12.2}, {632.4, -12.2}, {632.4, -32.2}, {637.4, -32.2}}, color = {0, 0, 127}));
      connect(powerCalc_pcc.P, filter3.u) annotation(
        Line(points = {{629.4, -1.56}, {637.4, -1.56}, {637.4, -0.56}}, color = {0, 0, 127}));
      connect(powerCalc_pcc.P, P_pupc.u1) annotation(
        Line(points = {{629.4, -1.56}, {631.4, -1.56}, {631.4, 54.44}, {680.4, 54.44}}, color = {0, 0, 127}));
      connect(realExpression43.y, clark3p4.V_abc_D[1]) annotation(
        Line(points = {{522, 279}, {541, 279}, {541, 274}, {548, 274}}, color = {0, 0, 127}));
      connect(ipc.q, ipq.u1) annotation(
        Line(points = {{622, 270}, {648, 270}, {648, 221.6}, {661, 221.6}}, color = {0, 0, 127}));
      connect(constant15.y, ipq.u2) annotation(
        Line(points = {{641, 217}, {661, 217}, {661, 214}}, color = {0, 0, 127}));
      connect(realExpression42.y, clark3p4.V_abc_D[2]) annotation(
        Line(points = {{522, 262}, {538, 262}, {538, 274}, {548, 274}}, color = {0, 0, 127}));
      connect(clark3p4.V_ab_D[1], ipc.alpha) annotation(
        Line(points = {{571, 274}, {581.5, 274}, {581.5, 280}, {599, 280}}, color = {0, 0, 127}));
      connect(realExpression41.y, clark3p4.V_abc_D[3]) annotation(
        Line(points = {{527, 243}, {539, 243}, {539, 274}, {548, 274}}, color = {0, 0, 127}));
      connect(ipc.d, ipd.u1) annotation(
        Line(points = {{622, 276}, {652, 276}, {652, 248}, {667, 248}}, color = {0, 0, 127}));
      connect(constant15.y, ipd.u2) annotation(
        Line(points = {{641, 217}, {654, 217}, {654, 238}, {667, 238}}, color = {0, 0, 127}));
      connect(ipc.beta, clark3p4.V_ab_D[2]) annotation(
        Line(points = {{598.8, 271.4}, {585.6, 271.4}, {585.6, 274}, {571, 274}}, color = {0, 0, 127}));
  connect(inversePark_fz.alpha, clarkInv.A) annotation(
        Line(points = {{485, 145}, {520, 145}, {520, 144}}, color = {0, 0, 127}));
  connect(inversePark_fz.beta, clarkInv.B) annotation(
        Line(points = {{485, 137}, {502.5, 137}, {502.5, 138}, {520, 138}}, color = {0, 0, 127}));
  connect(constant2.y, clarkInv.C) annotation(
        Line(points = {{500, 98}, {510, 98}, {510, 133}, {520, 133}}, color = {0, 0, 127}));
      annotation(
        Diagram(coordinateSystem(extent = {{-150, 300}, {700, -130}}, grid = {1, 1})),
        Icon(coordinateSystem(extent = {{-120, -120}, {480, 160}}, grid = {1, 1})));
    end ART22;

    model ART32
      parameter Real Rf = 3.53e-4;
      parameter Real Lf = 2.81e-5;
      parameter Real Rr = 1.47;
      parameter Real Lr = 7.32e-3;
      VSC_FZ.Testes.VSC vsc annotation(
        Placement(visible = true, transformation(origin = {-1, 27}, extent = {{-13, -13}, {13, 13}}, rotation = 0)));
      Modelica.Electrical.Analog.Sources.ConstantVoltage constantVoltage(V = 395.18) annotation(
        Placement(visible = true, transformation(origin = {-86, 18}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      Modelica.Electrical.Analog.Basic.Resistor resistor(R = Rf, alpha = 0) annotation(
        Placement(visible = true, transformation(origin = {87, 34}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Resistor resistor1(R = Rf, alpha = 0) annotation(
        Placement(visible = true, transformation(origin = {90, -26}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Resistor resistor2(R = Rf, alpha = 0) annotation(
        Placement(visible = true, transformation(origin = {92, -86}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Inductor inductor(L = Lf) annotation(
        Placement(visible = true, transformation(origin = {115, 34}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Inductor inductor1(L = Lf) annotation(
        Placement(visible = true, transformation(origin = {118, -26}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      VSC_FZ.Testes.Clark3p clark3p annotation(
        Placement(visible = true, transformation(origin = {-87, 129}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.VoltageSensor vsa annotation(
        Placement(visible = true, transformation(origin = {142, 53}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.VoltageSensor vsb annotation(
        Placement(visible = true, transformation(origin = {141, -6}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.VoltageSensor vscv annotation(
        Placement(visible = true, transformation(origin = {138, -64}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.CurrentSensor ia annotation(
        Placement(visible = true, transformation(origin = {54, 34}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.CurrentSensor ib annotation(
        Placement(visible = true, transformation(origin = {56, -26}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.CurrentSensor ic annotation(
        Placement(visible = true, transformation(origin = {58, -86}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
      VSC_FZ.Testes.Clark3p clark3p1 annotation(
        Placement(visible = true, transformation(origin = {205, 144}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression(y = ia.i) annotation(
        Placement(visible = true, transformation(origin = {171, 144}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression1(y = ib.i) annotation(
        Placement(visible = true, transformation(origin = {171, 128}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression2(y = ic.i) annotation(
        Placement(visible = true, transformation(origin = {171, 112}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression3(y = pLL_fz2.theta) annotation(
        Placement(visible = true, transformation(origin = {215.5, 116}, extent = {{-20.5, -9}, {20.5, 9}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression4(y = vsb.v) annotation(
        Placement(visible = true, transformation(origin = {-129, 123}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression5(y = vscv.v) annotation(
        Placement(visible = true, transformation(origin = {-129, 107}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression6(y = vsa.v) annotation(
        Placement(visible = true, transformation(origin = {-129, 139}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Sources.SineVoltage sineVoltage(V = 11.268e3, freqHz = 60) annotation(
        Placement(visible = true, transformation(origin = {362, 39}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Inductor inductor2(L = Lf) annotation(
        Placement(visible = true, transformation(origin = {120, -86}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Sources.SineVoltage sineVoltage1(V = 11.268e3, freqHz = 60, phase = -2.094395102393195) annotation(
        Placement(visible = true, transformation(origin = {364, -21}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Sources.SineVoltage sineVoltage2(V = 11.268e3, freqHz = 60, phase = 2.094395102393195) annotation(
        Placement(visible = true, transformation(origin = {364, -81}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Ground ground1 annotation(
        Placement(visible = true, transformation(origin = {399, -21}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
      Modelica.Blocks.Sources.RealExpression realExpression7(y = pLL_fz2.d) annotation(
        Placement(visible = true, transformation(origin = {344, 173.5}, extent = {{-17, -9.5}, {17, 9.5}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression8(y = pLL_fz2.q) annotation(
        Placement(visible = true, transformation(origin = {311.5, 98.5}, extent = {{-15.5, -9.5}, {15.5, 9.5}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression9(y = Controle_potencia.idref) annotation(
        Placement(visible = true, transformation(origin = {333, 158.5}, extent = {{-28, -9.5}, {28, 9.5}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression10(y = Correntes.d) annotation(
        Placement(visible = true, transformation(origin = {341, 145}, extent = {{-19, -10}, {19, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression11(y = Correntes.q) annotation(
        Placement(visible = true, transformation(origin = {341.5, 126}, extent = {{-19.5, -10}, {19.5, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression12(y = Controle_potencia.iqref) annotation(
        Placement(visible = true, transformation(origin = {335, 113.5}, extent = {{-26, -8.5}, {26, 8.5}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression13(y = pLL_fz2.theta) annotation(
        Placement(visible = true, transformation(origin = {445, 97.5}, extent = {{-23, -11.5}, {23, 11.5}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant2(k = 0) annotation(
        Placement(visible = true, transformation(origin = {489, 98}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression14(y = clarkInv.a) annotation(
        Placement(visible = true, transformation(origin = {-73, -56.5}, extent = {{-17, -8.5}, {17, 8.5}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression15(y = clarkInv.b) annotation(
        Placement(visible = true, transformation(origin = {-75, -83.5}, extent = {{-17, -9.5}, {17, 9.5}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression16(y = clarkInv.c) annotation(
        Placement(visible = true, transformation(origin = {-76, -111.5}, extent = {{-18, -9.5}, {18, 9.5}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.VoltageSensor vdc_2 annotation(
        Placement(visible = true, transformation(origin = {-109, 21}, extent = {{-7, -7}, {7, 7}}, rotation = -90)));
      Modelica.Blocks.Math.Division division annotation(
        Placement(visible = true, transformation(origin = {434, 157}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Math.Division division1 annotation(
        Placement(visible = true, transformation(origin = {432, 121}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Nonlinear.Limiter limiter(limitsAtInit = true, uMax = 1, uMin = -1) annotation(
        Placement(visible = true, transformation(origin = {-28, -58}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
      Modelica.Blocks.Nonlinear.Limiter limiter1(limitsAtInit = true, uMax = 1, uMin = -1) annotation(
        Placement(visible = true, transformation(origin = {-28, -84}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
      Modelica.Blocks.Nonlinear.Limiter limiter2(limitsAtInit = true, uMax = 1, uMin = -1) annotation(
        Placement(visible = true, transformation(origin = {-28, -112}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
      VSC_FZ.Testes.InnerControl2 Controle_corrente(Ictrl_KI = 0.5545, Ictrl_KP = 4.42e-2, lw = 2*3.1416*60*Lf, r = Rf) annotation(
        Placement(visible = true, transformation(origin = {394.5, 144.5}, extent = {{-14.5, -14.5}, {14.5, 14.5}}, rotation = 0)));
      Modelica.Blocks.Sources.Step Pref(height = 841.8e3, startTime = 0.1) annotation(
        Placement(visible = true, transformation(origin = {39, 139}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Step Qref(height = 358.6e3, startTime = 0.3) annotation(
        Placement(visible = true, transformation(origin = {37, 75}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      VSC_FZ.Testes.InversePark_fz inversePark_fz annotation(
        Placement(visible = true, transformation(origin = {474, 141}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      VSC_FZ.Testes.Park_Fz Correntes annotation(
        Placement(visible = true, transformation(origin = {243, 140}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      VSC_FZ.Testes.PLL_fz2 pLL_fz2(frequency = 60) annotation(
        Placement(visible = true, transformation(origin = {-34, 133}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Sources.ConstantVoltage constantVoltage1(V = 395.18) annotation(
        Placement(visible = true, transformation(origin = {-86, -16}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      Modelica.Electrical.Analog.Basic.Ground ground annotation(
        Placement(visible = true, transformation(origin = {-65, 70}, extent = {{-54, -86}, {-38, -70}}, rotation = 0)));
      Modelica.Blocks.Continuous.Filter filter(analogFilter = Modelica.Blocks.Types.AnalogFilter.Bessel, f_cut = 60, filterType = Modelica.Blocks.Types.FilterType.LowPass) annotation(
        Placement(visible = true, transformation(origin = {406, -18}, extent = {{84, 4}, {104, 24}}, rotation = 0)));
      VSC_FZ.Testes.PowerCalc powerCalc annotation(
        Placement(visible = true, transformation(origin = {465, -11}, extent = {{-14, -14}, {14, 14}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression17(y = vsa.v) annotation(
        Placement(visible = true, transformation(origin = {431, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression24(y = vsb.v) annotation(
        Placement(visible = true, transformation(origin = {431, -3}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression25(y = vscv.v) annotation(
        Placement(visible = true, transformation(origin = {431, -16}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression26(y = ia.i) annotation(
        Placement(visible = true, transformation(origin = {431, -31}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression27(y = ib.i) annotation(
        Placement(visible = true, transformation(origin = {431, -47}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression28(y = ic.i) annotation(
        Placement(visible = true, transformation(origin = {431, -63}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Continuous.Filter filter1(analogFilter = Modelica.Blocks.Types.AnalogFilter.Bessel, f_cut = 60, filterType = Modelica.Blocks.Types.FilterType.LowPass) annotation(
        Placement(visible = true, transformation(origin = {406, -49}, extent = {{84, 4}, {104, 24}}, rotation = 0)));
      VSC_FZ.Testes.PowerControl1 Controle_potencia(Pctrl_KI = 0.168, Pctrl_KP = 5.36e-4) annotation(
        Placement(visible = true, transformation(origin = {108.5, 104.5}, extent = {{-17.5, -17.5}, {17.5, 17.5}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression23(y = powerCalc.P) annotation(
        Placement(visible = true, transformation(origin = {33.5, 111}, extent = {{-20.5, -10}, {20.5, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression29(y = powerCalc.Q) annotation(
        Placement(visible = true, transformation(origin = {33.5, 94}, extent = {{-20.5, -10}, {20.5, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression18(y = vdc_2.v) annotation(
        Placement(visible = true, transformation(origin = {388, 78.5}, extent = {{-23, -11.5}, {23, 11.5}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Resistor resistor3(R = Rr, alpha = 0) annotation(
        Placement(visible = true, transformation(origin = {299, 39}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Resistor resistor4(R = Rr, alpha = 0) annotation(
        Placement(visible = true, transformation(origin = {298, -21}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Resistor resistor5(R = Rr, alpha = 0) annotation(
        Placement(visible = true, transformation(origin = {298, -81}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Ground ground2 annotation(
        Placement(visible = true, transformation(origin = {161.5, -6.5}, extent = {{-6.5, -6.5}, {6.5, 6.5}}, rotation = 90)));
      Modelica.Electrical.Analog.Basic.Ground ground3 annotation(
        Placement(visible = true, transformation(origin = {160.5, 52.5}, extent = {{-6.5, -6.5}, {6.5, 6.5}}, rotation = 90)));
      Modelica.Electrical.Analog.Basic.Ground ground4 annotation(
        Placement(visible = true, transformation(origin = {160.5, -64.5}, extent = {{-6.5, -6.5}, {6.5, 6.5}}, rotation = 90)));
      Modelica.Electrical.Analog.Basic.Inductor inductor3(L = Lr) annotation(
        Placement(visible = true, transformation(origin = {331, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Inductor inductor4(L = Lr) annotation(
        Placement(visible = true, transformation(origin = {330, -21}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Inductor inductor5(L = Lr) annotation(
        Placement(visible = true, transformation(origin = {330, -81}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      VSC_FZ.Testes.y_D_Transformer_ideal y_D_Transformer_ideal(Vac_primary = 440, Vac_secondary = 13800) annotation(
        Placement(visible = true, transformation(origin = {208, -28}, extent = {{-22, -22}, {22, 22}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Ground ground5 annotation(
        Placement(visible = true, transformation(origin = {302.5, 14.5}, extent = {{-6.5, -6.5}, {6.5, 6.5}}, rotation = 90)));
      Modelica.Electrical.Analog.Sensors.VoltageSensor pca annotation(
        Placement(visible = true, transformation(origin = {284, 15}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.VoltageSensor pcc annotation(
        Placement(visible = true, transformation(origin = {282, -106}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Ground ground6 annotation(
        Placement(visible = true, transformation(origin = {300.5, -106.5}, extent = {{-6.5, -6.5}, {6.5, 6.5}}, rotation = 90)));
      Modelica.Electrical.Analog.Basic.Ground ground7 annotation(
        Placement(visible = true, transformation(origin = {305.5, -44.5}, extent = {{-6.5, -6.5}, {6.5, 6.5}}, rotation = 90)));
      Modelica.Electrical.Analog.Sensors.VoltageSensor pcb annotation(
        Placement(visible = true, transformation(origin = {287, -44}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.CurrentSensor currentSensor annotation(
        Placement(visible = true, transformation(origin = {-97, 1}, extent = {{-4, -4}, {4, 4}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Ground ground8 annotation(
        Placement(visible = true, transformation(origin = {92.5, 9.5}, extent = {{-6.5, -6.5}, {6.5, 6.5}}, rotation = 90)));
      Modelica.Electrical.Analog.Sensors.VoltageSensor vta annotation(
        Placement(visible = true, transformation(origin = {74, 10}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.VoltageSensor vtb annotation(
        Placement(visible = true, transformation(origin = {78, -51}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Ground ground9 annotation(
        Placement(visible = true, transformation(origin = {96.5, -51.5}, extent = {{-6.5, -6.5}, {6.5, 6.5}}, rotation = 90)));
      Modelica.Electrical.Analog.Basic.Ground ground10 annotation(
        Placement(visible = true, transformation(origin = {101.5, -115.5}, extent = {{-6.5, -6.5}, {6.5, 6.5}}, rotation = 90)));
      Modelica.Electrical.Analog.Sensors.VoltageSensor vtc annotation(
        Placement(visible = true, transformation(origin = {83, -115}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Math.Division pcd annotation(
        Placement(visible = true, transformation(origin = {213, 210}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression30(y = pcc.v) annotation(
        Placement(visible = true, transformation(origin = {77, 190}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression31(y = pcb.v) annotation(
        Placement(visible = true, transformation(origin = {77, 206}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression32(y = pca.v) annotation(
        Placement(visible = true, transformation(origin = {77, 222}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Division pcq annotation(
        Placement(visible = true, transformation(origin = {212, 186}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant3(k = 11.267e3) annotation(
        Placement(visible = true, transformation(origin = {173, 185}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      VSC_FZ.Testes.Clark3p clark3p3 annotation(
        Placement(visible = true, transformation(origin = {121, 209}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      VSC_FZ.Testes.PLL_fz2 pLL_fz21(frequency = 60) annotation(
        Placement(visible = true, transformation(origin = {174, 213}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Division vsq annotation(
        Placement(visible = true, transformation(origin = {-16, 81}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Math.Division vsd annotation(
        Placement(visible = true, transformation(origin = {-15, 105}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant4(k = 359.258) annotation(
        Placement(visible = true, transformation(origin = {-55, 80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Division Pref_pu annotation(
        Placement(visible = true, transformation(origin = {130, 159}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant5(k = 915e3) annotation(
        Placement(visible = true, transformation(origin = {97, 143}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant6(k = 915e3) annotation(
        Placement(visible = true, transformation(origin = {506, 31}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Division P_pu annotation(
        Placement(visible = true, transformation(origin = {539, 47}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant7(k = 915e3) annotation(
        Placement(visible = true, transformation(origin = {506, -88}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Division Q_pu annotation(
        Placement(visible = true, transformation(origin = {539, -72}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant8(k = 1698) annotation(
        Placement(visible = true, transformation(origin = {263, 173}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Division id_pu annotation(
        Placement(visible = true, transformation(origin = {296, 189}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant9(k = 1698) annotation(
        Placement(visible = true, transformation(origin = {238, 95}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Division iq_pu annotation(
        Placement(visible = true, transformation(origin = {271, 111}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant10(k = 1698) annotation(
        Placement(visible = true, transformation(origin = {356, 205}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Division idref_pu annotation(
        Placement(visible = true, transformation(origin = {389, 221}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression19(y = Controle_potencia.idref) annotation(
        Placement(visible = true, transformation(origin = {326, 232.5}, extent = {{-28, -9.5}, {28, 9.5}}, rotation = 0)));
      Modelica.Blocks.Math.Division iqref_pu annotation(
        Placement(visible = true, transformation(origin = {503, 221}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant11(k = 1698) annotation(
        Placement(visible = true, transformation(origin = {470, 205}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression21(y = Controle_potencia.iqref) annotation(
        Placement(visible = true, transformation(origin = {444, 229.5}, extent = {{-26, -8.5}, {26, 8.5}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant12(k = 915e3) annotation(
        Placement(visible = true, transformation(origin = {116, 259}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Division Qref_pu annotation(
        Placement(visible = true, transformation(origin = {149, 275}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression20(y = Qref.y) annotation(
        Placement(visible = true, transformation(origin = {107, 283.5}, extent = {{-17, -9.5}, {17, 9.5}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.CurrentSensor ipca annotation(
        Placement(visible = true, transformation(origin = {260, 39}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.CurrentSensor ipcb annotation(
        Placement(visible = true, transformation(origin = {262, -21}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.CurrentSensor ipcc annotation(
        Placement(visible = true, transformation(origin = {264, -81}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression39(y = ia.i) annotation(
        Placement(visible = true, transformation(origin = {-118, 245}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression40(y = ib.i) annotation(
        Placement(visible = true, transformation(origin = {-118, 229}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      VSC_FZ.Testes.PLL_fz2 correntes_tot(frequency = 60) annotation(
        Placement(visible = true, transformation(origin = {-18, 243}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant14(k = 1698) annotation(
        Placement(visible = true, transformation(origin = {-39, 190}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Division itot_d annotation(
        Placement(visible = true, transformation(origin = {1, 215}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Math.Division itot_q annotation(
        Placement(visible = true, transformation(origin = {0, 191}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      VSC_FZ.Testes.Clark3p clark3p2 annotation(
        Placement(visible = true, transformation(origin = {-71, 239}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression38(y = ic.i) annotation(
        Placement(visible = true, transformation(origin = {-118, 213}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant1(k = 915e3) annotation(
        Placement(visible = true, transformation(origin = {655, 34}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression22(y = ipcb.i) annotation(
        Placement(visible = true, transformation(origin = {580, -44}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant13(k = 915e3) annotation(
        Placement(visible = true, transformation(origin = {655, -85}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Continuous.Filter filter2(analogFilter = Modelica.Blocks.Types.AnalogFilter.Bessel, f_cut = 60, filterType = Modelica.Blocks.Types.FilterType.LowPass) annotation(
        Placement(visible = true, transformation(origin = {555, -46}, extent = {{84, 4}, {104, 24}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression33(y = ipca.i) annotation(
        Placement(visible = true, transformation(origin = {580, -28}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Division P_pupc annotation(
        Placement(visible = true, transformation(origin = {688, 50}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Continuous.Filter filter3(analogFilter = Modelica.Blocks.Types.AnalogFilter.Bessel, f_cut = 60, filterType = Modelica.Blocks.Types.FilterType.LowPass) annotation(
        Placement(visible = true, transformation(origin = {555, -15}, extent = {{84, 4}, {104, 24}}, rotation = 0)));
      Modelica.Blocks.Math.Division Q_pupc annotation(
        Placement(visible = true, transformation(origin = {688, -69}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression34(y = pcb.v) annotation(
        Placement(visible = true, transformation(origin = {580, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression35(y = ipcc.i) annotation(
        Placement(visible = true, transformation(origin = {580, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression36(y = pcc.v) annotation(
        Placement(visible = true, transformation(origin = {580, -13}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      VSC_FZ.Testes.PowerCalc powerCalc_pcc annotation(
        Placement(visible = true, transformation(origin = {614, -8}, extent = {{-14, -14}, {14, 14}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression37(y = pca.v) annotation(
        Placement(visible = true, transformation(origin = {580, 13}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      VSC_FZ.Testes.PLL_fz2 ipc(frequency = 60) annotation(
        Placement(visible = true, transformation(origin = {611, 276}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression41(y = ipcc.i) annotation(
        Placement(visible = true, transformation(origin = {516, 243}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Division ipd annotation(
        Placement(visible = true, transformation(origin = {677, 243}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression42(y = ipcb.i) annotation(
        Placement(visible = true, transformation(origin = {511, 262}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      VSC_FZ.Testes.Clark3p clark3p4 annotation(
        Placement(visible = true, transformation(origin = {560, 274}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant15(k = 54.13) annotation(
        Placement(visible = true, transformation(origin = {630, 217}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression43(y = ipca.i) annotation(
        Placement(visible = true, transformation(origin = {511, 279}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Division ipq annotation(
        Placement(visible = true, transformation(origin = {669, 218}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression44(y = ib.i) annotation(
        Placement(visible = true, transformation(origin = {220, 264}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression45(y = ic.i) annotation(
        Placement(visible = true, transformation(origin = {220, 248}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression46(y = ia.i) annotation(
        Placement(visible = true, transformation(origin = {220, 280}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Division ia_pu annotation(
        Placement(visible = true, transformation(origin = {277, 289}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Math.Division ib_pu annotation(
        Placement(visible = true, transformation(origin = {277, 262}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Math.Division ic_pu annotation(
        Placement(visible = true, transformation(origin = {277, 236}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant16(k = 1698) annotation(
        Placement(visible = true, transformation(origin = {240, 217}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  ClarkInv clarkInv annotation(
        Placement(visible = true, transformation(origin = {523, 138}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(constantVoltage.p, vsc.pin_p) annotation(
        Line(points = {{-86, 28}, {-86, 34}, {-15, 34}, {-15, 33.5}}, color = {0, 0, 255}));
      connect(resistor.n, inductor.p) annotation(
        Line(points = {{97, 34}, {105, 34}}, color = {0, 0, 255}));
      connect(resistor1.n, inductor1.p) annotation(
        Line(points = {{100, -26}, {108, -26}}, color = {0, 0, 255}));
      connect(vsc.Vta, ia.p) annotation(
        Line(points = {{13, 35}, {30, 35}, {30, 34}, {48, 34}}, color = {0, 0, 255}));
      connect(ia.n, resistor.p) annotation(
        Line(points = {{60, 34}, {77, 34}}, color = {0, 0, 255}));
      connect(vsc.Vtb, ib.p) annotation(
        Line(points = {{13, 27}, {40, 27}, {40, -26}, {50, -26}}, color = {0, 0, 255}));
      connect(ib.n, resistor1.p) annotation(
        Line(points = {{62, -26}, {80, -26}}, color = {0, 0, 255}));
      connect(ic.n, resistor2.p) annotation(
        Line(points = {{64, -86}, {82, -86}}, color = {0, 0, 255}));
      connect(ic.p, vsc.Vtc) annotation(
        Line(points = {{52, -86}, {24, -86}, {24, 19}, {13, 19}}, color = {0, 0, 255}));
      connect(realExpression.y, clark3p1.V_abc_D[1]) annotation(
        Line(points = {{182, 144}, {192, 144}}, color = {0, 0, 127}));
      connect(realExpression1.y, clark3p1.V_abc_D[2]) annotation(
        Line(points = {{182, 128}, {184, 128}, {184, 144}, {192, 144}}, color = {0, 0, 127}));
      connect(realExpression2.y, clark3p1.V_abc_D[3]) annotation(
        Line(points = {{182, 112}, {184, 112}, {184, 144}, {192, 144}}, color = {0, 0, 127}));
      connect(realExpression6.y, clark3p.V_abc_D[1]) annotation(
        Line(points = {{-118, 139}, {-104, 139}, {-104, 129}, {-100, 129}}, color = {0, 0, 127}));
      connect(realExpression4.y, clark3p.V_abc_D[2]) annotation(
        Line(points = {{-118, 123}, {-112, 123}, {-112, 129}, {-100, 129}}, color = {0, 0, 127}));
      connect(realExpression5.y, clark3p.V_abc_D[3]) annotation(
        Line(points = {{-118, 107}, {-112, 107}, {-112, 129}, {-100, 129}}, color = {0, 0, 127}));
      connect(inductor2.p, resistor2.n) annotation(
        Line(points = {{110, -86}, {102, -86}}, color = {0, 0, 255}));
      connect(sineVoltage.n, ground1.p) annotation(
        Line(points = {{372, 39}, {389, 39}, {389, -21}}, color = {0, 0, 255}));
      connect(sineVoltage1.n, ground1.p) annotation(
        Line(points = {{374, -21}, {389, -21}}, color = {0, 0, 255}));
      connect(sineVoltage2.n, ground1.p) annotation(
        Line(points = {{374, -81}, {389, -81}, {389, -21}}, color = {0, 0, 255}));
      connect(vsa.p, inductor.n) annotation(
        Line(points = {{135, 53}, {125, 53}, {125, 34}}, color = {0, 0, 255}));
      connect(vsb.p, inductor1.n) annotation(
        Line(points = {{134, -6}, {128, -6}, {128, -26}}, color = {0, 0, 255}));
      connect(vscv.p, inductor2.n) annotation(
        Line(points = {{130, -64}, {130, -86}}, color = {0, 0, 255}));
      connect(vdc_2.p, constantVoltage.p) annotation(
        Line(points = {{-109, 28}, {-87, 28}}, color = {0, 0, 255}));
      connect(vdc_2.n, constantVoltage.n) annotation(
        Line(points = {{-109, 14}, {-109, 8}, {-87, 8}}, color = {0, 0, 255}));
      connect(realExpression14.y, limiter.u) annotation(
        Line(points = {{-54, -56.5}, {-46, -56.5}, {-46, -58}, {-38, -58}}, color = {0, 0, 127}));
      connect(limiter.y, vsc.ma) annotation(
        Line(points = {{-19, -58}, {-8, -58}, {-8, 12}}, color = {0, 0, 127}));
      connect(realExpression15.y, limiter1.u) annotation(
        Line(points = {{-56, -83.5}, {-47, -83.5}, {-47, -84}, {-38, -84}}, color = {0, 0, 127}));
      connect(realExpression16.y, limiter2.u) annotation(
        Line(points = {{-56, -111.5}, {-47, -111.5}, {-47, -112}, {-38, -112}}, color = {0, 0, 127}));
      connect(limiter1.y, vsc.mb) annotation(
        Line(points = {{-20, -84}, {0, -84}, {0, 12}}, color = {0, 0, 127}));
      connect(limiter2.y, vsc.mc) annotation(
        Line(points = {{-20, -112}, {6, -112}, {6, 12}}, color = {0, 0, 127}));
      connect(Controle_corrente.vdref, division.u1) annotation(
        Line(points = {{410, 153}, {410, 161.22}, {424.45, 161.22}}, color = {0, 0, 127}));
      connect(Controle_corrente.vqref, division1.u1) annotation(
        Line(points = {{410, 138}, {410, 125.23}, {422.45, 125.23}}, color = {0, 0, 127}));
      connect(realExpression7.y, Controle_corrente.vd) annotation(
        Line(points = {{362.7, 173.5}, {373.7, 173.5}, {373.7, 158}, {379, 158}}, color = {0, 0, 127}));
      connect(realExpression10.y, Controle_corrente.Id) annotation(
        Line(points = {{361.9, 145}, {370.4, 145}, {370.4, 146}, {379, 146}}, color = {0, 0, 127}));
      connect(realExpression11.y, Controle_corrente.Iq) annotation(
        Line(points = {{362.95, 126}, {369.95, 126}, {369.95, 141}, {379, 141}}, color = {0, 0, 127}));
      connect(realExpression8.y, Controle_corrente.vq) annotation(
        Line(points = {{328.55, 98.5}, {379, 98.5}, {379, 131}}, color = {0, 0, 127}));
      connect(division.y, inversePark_fz.d) annotation(
        Line(points = {{441.7, 157}, {453.7, 157}, {453.7, 145}, {461.7, 145}}, color = {0, 0, 127}));
      connect(division1.y, inversePark_fz.q) annotation(
        Line(points = {{439.7, 121}, {455.7, 121}, {455.7, 137}, {461.7, 137}}, color = {0, 0, 127}));
      connect(inversePark_fz.theta, realExpression13.y) annotation(
        Line(points = {{474, 129}, {474, 97.5}, {470, 97.5}}, color = {0, 0, 127}));
      connect(clark3p1.V_ab_D[1], Correntes.alpha) annotation(
        Line(points = {{216, 144}, {231, 144}}, color = {0, 0, 127}));
      connect(clark3p1.V_ab_D[2], Correntes.beta) annotation(
        Line(points = {{216, 144}, {224, 144}, {224, 136}, {231, 136}}, color = {0, 0, 127}));
      connect(Correntes.theta, realExpression3.y) annotation(
        Line(points = {{243, 128}, {243, 116}, {238, 116}}, color = {0, 0, 127}));
      connect(clark3p.V_ab_D[1], pLL_fz2.alpha) annotation(
        Line(points = {{-76, 129}, {-63.5, 129}, {-63.5, 137}, {-46, 137}}, color = {0, 0, 127}));
      connect(pLL_fz2.beta, clark3p.V_ab_D[2]) annotation(
        Line(points = {{-46.2, 128.4}, {-59.2, 128.4}, {-59.2, 129.4}, {-76.2, 129.4}}, color = {0, 0, 127}));
      connect(constantVoltage1.n, vsc.pin_n) annotation(
        Line(points = {{-86, -26}, {-60, -26}, {-60, 20.5}, {-15, 20.5}}, color = {0, 0, 255}));
      connect(constantVoltage.n, constantVoltage1.p) annotation(
        Line(points = {{-86, 8}, {-86, -6}}, color = {0, 0, 255}));
      connect(realExpression17.y, powerCalc.V1) annotation(
        Line(points = {{442, 10}, {447, 10}, {447, 1}, {450, 1}}, color = {0, 0, 127}));
      connect(realExpression24.y, powerCalc.V2) annotation(
        Line(points = {{442, -3}, {450, -3}}, color = {0, 0, 127}));
      connect(realExpression25.y, powerCalc.V3) annotation(
        Line(points = {{442, -16}, {444, -16}, {444, -7}, {450, -7}}, color = {0, 0, 127}));
      connect(realExpression26.y, powerCalc.I1) annotation(
        Line(points = {{442, -31}, {446, -31}, {446, -14}, {450, -14}}, color = {0, 0, 127}));
      connect(realExpression27.y, powerCalc.I2) annotation(
        Line(points = {{442, -47}, {447, -47}, {447, -18}, {450, -18}}, color = {0, 0, 127}));
      connect(realExpression28.y, powerCalc.I3) annotation(
        Line(points = {{442, -63}, {448, -63}, {448, -23}, {450, -23}}, color = {0, 0, 127}));
      connect(powerCalc.P, filter.u) annotation(
        Line(points = {{480.4, -4.56}, {488.4, -4.56}, {488.4, -3.56}}, color = {0, 0, 127}));
      connect(powerCalc.Q, filter1.u) annotation(
        Line(points = {{480.4, -15.2}, {483.4, -15.2}, {483.4, -35.2}, {488.4, -35.2}}, color = {0, 0, 127}));
      connect(Qref.y, Controle_potencia.Q_ref) annotation(
        Line(points = {{48, 75}, {78, 75}, {78, 94}, {89, 94}}, color = {0, 0, 127}));
      connect(Pref.y, Controle_potencia.P_ref) annotation(
        Line(points = {{50, 139}, {81, 139}, {81, 113}, {89, 113}}, color = {0, 0, 127}));
      connect(realExpression9.y, Controle_corrente.Id_ref) annotation(
        Line(points = {{364, 158.5}, {369.25, 158.5}, {369.25, 152}, {379, 152}}, color = {0, 0, 127}));
      connect(realExpression12.y, Controle_corrente.Iqref) annotation(
        Line(points = {{364, 113.5}, {374.1, 113.5}, {374.1, 136}, {379, 136}}, color = {0, 0, 127}));
      connect(realExpression29.y, Controle_potencia.Q) annotation(
        Line(points = {{56, 94}, {72, 94}, {72, 100}, {89, 100}}, color = {0, 0, 127}));
      connect(realExpression23.y, Controle_potencia.P) annotation(
        Line(points = {{56, 111}, {67, 111}, {67, 106}, {89, 106}}, color = {0, 0, 127}));
      connect(realExpression18.y, division1.u2) annotation(
        Line(points = {{413.3, 78.5}, {413.3, 117}, {424.3, 117}}, color = {0, 0, 127}));
      connect(realExpression18.y, division.u2) annotation(
        Line(points = {{413.3, 78.5}, {413.3, 153}, {426.3, 153}}, color = {0, 0, 127}));
      connect(vsb.n, ground2.p) annotation(
        Line(points = {{148, -6}, {155, -6}}, color = {0, 0, 255}));
      connect(vsa.n, ground3.p) annotation(
        Line(points = {{149, 53}, {154, 53}}, color = {0, 0, 255}));
      connect(vscv.n, ground4.p) annotation(
        Line(points = {{146, -64}, {154, -64}}, color = {0, 0, 255}));
      connect(resistor4.n, inductor4.n) annotation(
        Line(points = {{308, -21}, {320, -21}}, color = {0, 0, 255}));
      connect(resistor5.n, inductor5.n) annotation(
        Line(points = {{308, -81}, {320, -81}}, color = {0, 0, 255}));
      connect(inductor5.p, sineVoltage2.p) annotation(
        Line(points = {{340, -81}, {354, -81}}, color = {0, 0, 255}));
      connect(inductor4.p, sineVoltage1.p) annotation(
        Line(points = {{340, -21}, {354, -21}}, color = {0, 0, 255}));
      connect(resistor3.n, inductor3.p) annotation(
        Line(points = {{309, 39}, {313, 39}, {313, 40}, {321, 40}}, color = {0, 0, 255}));
      connect(inductor3.n, sineVoltage.p) annotation(
        Line(points = {{341, 40}, {352, 40}, {352, 39}}, color = {0, 0, 255}));
      connect(y_D_Transformer_ideal.pin_p, inductor.n) annotation(
        Line(points = {{184, -14}, {180, -14}, {180, 34}, {125, 34}}, color = {0, 0, 255}));
      connect(y_D_Transformer_ideal.pin_p1, inductor1.n) annotation(
        Line(points = {{184, -25}, {128, -25}, {128, -26}}, color = {0, 0, 255}));
      connect(y_D_Transformer_ideal.pin_p2, inductor2.n) annotation(
        Line(points = {{184, -35}, {171, -35}, {171, -86}, {130, -86}}, color = {0, 0, 255}));
      connect(pca.n, ground5.p) annotation(
        Line(points = {{291, 15}, {296, 15}}, color = {0, 0, 255}));
      connect(pcc.n, ground6.p) annotation(
        Line(points = {{289, -106}, {294, -106}}, color = {0, 0, 255}));
      connect(pcb.n, ground7.p) annotation(
        Line(points = {{294, -44}, {299, -44}}, color = {0, 0, 255}));
      connect(pca.p, y_D_Transformer_ideal.pin_n) annotation(
        Line(points = {{277, 15}, {238, 15}, {238, -15}, {232, -15}}, color = {0, 0, 255}));
      connect(pcb.p, y_D_Transformer_ideal.pin_n1) annotation(
        Line(points = {{280, -44}, {248, -44}, {248, -24}, {232, -24}}, color = {0, 0, 255}));
      connect(pcc.p, y_D_Transformer_ideal.pin_n2) annotation(
        Line(points = {{275, -106}, {243, -106}, {243, -34}, {232, -34}}, color = {0, 0, 255}));
      connect(ground.p, currentSensor.p) annotation(
        Line(points = {{-111, 0}, {-101, 0}, {-101, 1}}, color = {0, 0, 255}));
      connect(currentSensor.n, constantVoltage.n) annotation(
        Line(points = {{-93, 1}, {-86, 1}, {-86, 8}}, color = {0, 0, 255}));
      connect(vta.n, ground8.p) annotation(
        Line(points = {{81, 10}, {86, 10}}, color = {0, 0, 255}));
      connect(vtb.n, ground9.p) annotation(
        Line(points = {{85, -51}, {90, -51}}, color = {0, 0, 255}));
      connect(vtc.n, ground10.p) annotation(
        Line(points = {{90, -115}, {95, -115}}, color = {0, 0, 255}));
      connect(vtc.p, resistor2.p) annotation(
        Line(points = {{76, -115}, {70, -115}, {70, -86}, {82, -86}}, color = {0, 0, 255}));
      connect(vtb.p, resistor1.p) annotation(
        Line(points = {{71, -51}, {65, -51}, {65, -26}, {80, -26}}, color = {0, 0, 255}));
      connect(vta.p, resistor.p) annotation(
        Line(points = {{67, 10}, {64, 10}, {64, 34}, {77, 34}}, color = {0, 0, 255}));
      connect(realExpression31.y, clark3p3.V_abc_D[2]) annotation(
        Line(points = {{88, 206}, {94, 206}, {94, 209}, {109, 209}}, color = {0, 0, 127}));
      connect(pLL_fz21.q, pcq.u1) annotation(
        Line(points = {{185, 206.6}, {193, 206.6}, {193, 189.6}, {204, 189.6}}, color = {0, 0, 127}));
      connect(realExpression32.y, clark3p3.V_abc_D[1]) annotation(
        Line(points = {{88, 222}, {102, 222}, {102, 209}, {109, 209}}, color = {0, 0, 127}));
      connect(realExpression30.y, clark3p3.V_abc_D[3]) annotation(
        Line(points = {{88, 190}, {94, 190}, {94, 209}, {109, 209}}, color = {0, 0, 127}));
      connect(pLL_fz21.d, pcd.u1) annotation(
        Line(points = {{185, 213}, {205, 213}, {205, 214}}, color = {0, 0, 127}));
      connect(clark3p3.V_ab_D[1], pLL_fz21.alpha) annotation(
        Line(points = {{132, 209}, {144.5, 209}, {144.5, 217}, {162, 217}}, color = {0, 0, 127}));
      connect(constant3.y, pcq.u2) annotation(
        Line(points = {{184, 185}, {204, 185}, {204, 182}}, color = {0, 0, 127}));
      connect(pLL_fz21.beta, clark3p3.V_ab_D[2]) annotation(
        Line(points = {{161.8, 208.4}, {148.8, 208.4}, {148.8, 209.4}, {131.8, 209.4}}, color = {0, 0, 127}));
      connect(constant3.y, pcd.u2) annotation(
        Line(points = {{184, 185}, {197, 185}, {197, 206}, {205, 206}}, color = {0, 0, 127}));
      connect(constant4.y, vsq.u2) annotation(
        Line(points = {{-44, 80}, {-24, 80}, {-24, 77}}, color = {0, 0, 127}));
      connect(constant4.y, vsd.u2) annotation(
        Line(points = {{-44, 80}, {-31, 80}, {-31, 101}, {-23, 101}}, color = {0, 0, 127}));
      connect(pLL_fz2.d, vsd.u1) annotation(
        Line(points = {{-23, 133}, {-33, 133}, {-33, 109}, {-23, 109}}, color = {0, 0, 127}));
      connect(pLL_fz2.q, vsq.u1) annotation(
        Line(points = {{-23, 127}, {-37, 127}, {-37, 85}, {-24, 85}}, color = {0, 0, 127}));
      connect(Pref_pu.u1, Pref.y) annotation(
        Line(points = {{122, 163}, {81, 163}, {81, 139}, {50, 139}}, color = {0, 0, 127}));
      connect(constant5.y, Pref_pu.u2) annotation(
        Line(points = {{108, 143}, {113, 143}, {113, 155}, {122, 155}}, color = {0, 0, 127}));
      connect(constant6.y, P_pu.u2) annotation(
        Line(points = {{517, 31}, {522, 31}, {522, 43}, {531, 43}}, color = {0, 0, 127}));
      connect(powerCalc.P, P_pu.u1) annotation(
        Line(points = {{480, -5}, {482, -5}, {482, 51}, {531, 51}}, color = {0, 0, 127}));
      connect(constant7.y, Q_pu.u2) annotation(
        Line(points = {{517, -88}, {522, -88}, {522, -76}, {531, -76}}, color = {0, 0, 127}));
      connect(powerCalc.Q, Q_pu.u1) annotation(
        Line(points = {{480, -15}, {484, -15}, {484, -68}, {531, -68}}, color = {0, 0, 127}));
      connect(constant8.y, id_pu.u2) annotation(
        Line(points = {{274, 173}, {279, 173}, {279, 185}, {288, 185}}, color = {0, 0, 127}));
      connect(id_pu.u1, Correntes.d) annotation(
        Line(points = {{288, 193}, {242, 193}, {242, 152}, {263, 152}, {263, 144}, {254, 144}}, color = {0, 0, 127}));
      connect(constant9.y, iq_pu.u2) annotation(
        Line(points = {{249, 95}, {254, 95}, {254, 107}, {263, 107}}, color = {0, 0, 127}));
      connect(iq_pu.u1, Correntes.q) annotation(
        Line(points = {{263, 115}, {257, 115}, {257, 136}, {254, 136}}, color = {0, 0, 127}));
      connect(constant10.y, idref_pu.u2) annotation(
        Line(points = {{367, 205}, {372, 205}, {372, 217}, {381, 217}}, color = {0, 0, 127}));
      connect(realExpression19.y, idref_pu.u1) annotation(
        Line(points = {{357, 233}, {375, 233}, {375, 225}, {381, 225}}, color = {0, 0, 127}));
      connect(constant11.y, iqref_pu.u2) annotation(
        Line(points = {{481, 205}, {486, 205}, {486, 217}, {495, 217}}, color = {0, 0, 127}));
      connect(realExpression21.y, iqref_pu.u1) annotation(
        Line(points = {{473, 229.5}, {484, 229.5}, {484, 225}, {495, 225}}, color = {0, 0, 127}));
      connect(constant12.y, Qref_pu.u2) annotation(
        Line(points = {{127, 259}, {132, 259}, {132, 271}, {141, 271}}, color = {0, 0, 127}));
      connect(realExpression20.y, Qref_pu.u1) annotation(
        Line(points = {{126, 284}, {133, 284}, {133, 279}, {141, 279}}, color = {0, 0, 127}));
      connect(resistor3.p, ipca.n) annotation(
        Line(points = {{289, 39}, {266, 39}}, color = {0, 0, 255}));
      connect(ipca.p, y_D_Transformer_ideal.pin_n) annotation(
        Line(points = {{254, 39}, {238, 39}, {238, -15}, {232, -15}}, color = {0, 0, 255}));
      connect(resistor4.p, ipcb.n) annotation(
        Line(points = {{288, -21}, {268, -21}}, color = {0, 0, 255}));
      connect(ipcb.p, y_D_Transformer_ideal.pin_n1) annotation(
        Line(points = {{256, -21}, {248, -21}, {248, -24}, {232, -24}}, color = {0, 0, 255}));
      connect(resistor5.p, ipcc.n) annotation(
        Line(points = {{288, -81}, {270, -81}}, color = {0, 0, 255}));
      connect(ipcc.p, y_D_Transformer_ideal.pin_n2) annotation(
        Line(points = {{258, -81}, {243, -81}, {243, -34}, {232, -34}}, color = {0, 0, 255}));
      connect(constant14.y, itot_q.u2) annotation(
        Line(points = {{-28, 190}, {-8, 190}, {-8, 187}}, color = {0, 0, 127}));
      connect(realExpression39.y, clark3p2.V_abc_D[1]) annotation(
        Line(points = {{-107, 245}, {-91, 245}, {-91, 239}, {-83, 239}}, color = {0, 0, 127}));
      connect(clark3p2.V_ab_D[1], correntes_tot.alpha) annotation(
        Line(points = {{-60, 239}, {-47.5, 239}, {-47.5, 247}, {-30, 247}}, color = {0, 0, 127}));
      connect(correntes_tot.q, itot_q.u1) annotation(
        Line(points = {{-7, 236.6}, {-21, 236.6}, {-21, 195}, {-8, 195}}, color = {0, 0, 127}));
      connect(constant14.y, itot_d.u2) annotation(
        Line(points = {{-28, 190}, {-15, 190}, {-15, 211}, {-7, 211}}, color = {0, 0, 127}));
      connect(realExpression40.y, clark3p2.V_abc_D[2]) annotation(
        Line(points = {{-107, 229}, {-91, 229}, {-91, 239}, {-83, 239}}, color = {0, 0, 127}));
      connect(correntes_tot.beta, clark3p2.V_ab_D[2]) annotation(
        Line(points = {{-30.2, 238.4}, {-43.2, 238.4}, {-43.2, 239.4}, {-60.2, 239.4}}, color = {0, 0, 127}));
      connect(realExpression38.y, clark3p2.V_abc_D[3]) annotation(
        Line(points = {{-107, 213}, {-91, 213}, {-91, 239}, {-83, 239}}, color = {0, 0, 127}));
      connect(correntes_tot.d, itot_d.u1) annotation(
        Line(points = {{-7, 243}, {-17, 243}, {-17, 219}, {-7, 219}}, color = {0, 0, 127}));
      connect(constant1.y, P_pupc.u2) annotation(
        Line(points = {{666, 34}, {671, 34}, {671, 46}, {680, 46}}, color = {0, 0, 127}));
      connect(realExpression33.y, powerCalc_pcc.I1) annotation(
        Line(points = {{591, -28}, {595, -28}, {595, -11}, {599, -11}}, color = {0, 0, 127}));
      connect(realExpression35.y, powerCalc_pcc.I3) annotation(
        Line(points = {{591, -60}, {597, -60}, {597, -20}, {599, -20}}, color = {0, 0, 127}));
      connect(constant13.y, Q_pupc.u2) annotation(
        Line(points = {{666, -85}, {671, -85}, {671, -73}, {680, -73}}, color = {0, 0, 127}));
      connect(realExpression37.y, powerCalc_pcc.V1) annotation(
        Line(points = {{591, 13}, {596, 13}, {596, 4}, {599, 4}}, color = {0, 0, 127}));
      connect(powerCalc_pcc.Q, Q_pupc.u1) annotation(
        Line(points = {{629.4, -12.2}, {633.4, -12.2}, {633.4, -65.2}, {680.4, -65.2}}, color = {0, 0, 127}));
      connect(realExpression22.y, powerCalc_pcc.I2) annotation(
        Line(points = {{591, -44}, {596, -44}, {596, -15}, {599, -15}}, color = {0, 0, 127}));
      connect(realExpression34.y, powerCalc_pcc.V2) annotation(
        Line(points = {{591, 0}, {599, 0}}, color = {0, 0, 127}));
      connect(realExpression36.y, powerCalc_pcc.V3) annotation(
        Line(points = {{591, -13}, {593, -13}, {593, -4}, {599, -4}}, color = {0, 0, 127}));
      connect(powerCalc_pcc.Q, filter2.u) annotation(
        Line(points = {{629.4, -12.2}, {632.4, -12.2}, {632.4, -32.2}, {637.4, -32.2}}, color = {0, 0, 127}));
      connect(powerCalc_pcc.P, filter3.u) annotation(
        Line(points = {{629.4, -1.56}, {637.4, -1.56}, {637.4, -0.56}}, color = {0, 0, 127}));
      connect(powerCalc_pcc.P, P_pupc.u1) annotation(
        Line(points = {{629.4, -1.56}, {631.4, -1.56}, {631.4, 54.44}, {680.4, 54.44}}, color = {0, 0, 127}));
      connect(realExpression43.y, clark3p4.V_abc_D[1]) annotation(
        Line(points = {{522, 279}, {541, 279}, {541, 274}, {548, 274}}, color = {0, 0, 127}));
      connect(ipc.q, ipq.u1) annotation(
        Line(points = {{622, 270}, {648, 270}, {648, 221.6}, {661, 221.6}}, color = {0, 0, 127}));
      connect(constant15.y, ipq.u2) annotation(
        Line(points = {{641, 217}, {661, 217}, {661, 214}}, color = {0, 0, 127}));
      connect(realExpression42.y, clark3p4.V_abc_D[2]) annotation(
        Line(points = {{522, 262}, {538, 262}, {538, 274}, {548, 274}}, color = {0, 0, 127}));
      connect(clark3p4.V_ab_D[1], ipc.alpha) annotation(
        Line(points = {{571, 274}, {581.5, 274}, {581.5, 280}, {599, 280}}, color = {0, 0, 127}));
      connect(realExpression41.y, clark3p4.V_abc_D[3]) annotation(
        Line(points = {{527, 243}, {539, 243}, {539, 274}, {548, 274}}, color = {0, 0, 127}));
      connect(ipc.d, ipd.u1) annotation(
        Line(points = {{622, 276}, {652, 276}, {652, 248}, {667, 248}}, color = {0, 0, 127}));
      connect(constant15.y, ipd.u2) annotation(
        Line(points = {{641, 217}, {654, 217}, {654, 238}, {667, 238}}, color = {0, 0, 127}));
      connect(ipc.beta, clark3p4.V_ab_D[2]) annotation(
        Line(points = {{598.8, 271.4}, {585.6, 271.4}, {585.6, 274}, {571, 274}}, color = {0, 0, 127}));
      connect(constant16.y, ic_pu.u2) annotation(
        Line(points = {{251, 217}, {259, 217}, {259, 232}, {269, 232}}, color = {0, 0, 127}));
      connect(constant16.y, ib_pu.u2) annotation(
        Line(points = {{251, 217}, {259, 217}, {259, 258}, {269, 258}}, color = {0, 0, 127}));
      connect(constant16.y, ia_pu.u2) annotation(
        Line(points = {{251, 217}, {259, 217}, {259, 285}, {269, 285}}, color = {0, 0, 127}));
      connect(realExpression46.y, ia_pu.u1) annotation(
        Line(points = {{231, 280}, {246, 280}, {246, 293}, {269, 293}}, color = {0, 0, 127}));
      connect(realExpression44.y, ib_pu.u1) annotation(
        Line(points = {{231, 264}, {269, 264}, {269, 266}}, color = {0, 0, 127}));
      connect(realExpression45.y, ic_pu.u1) annotation(
        Line(points = {{231, 248}, {262, 248}, {262, 240}, {269, 240}}, color = {0, 0, 127}));
  connect(inversePark_fz.alpha, clarkInv.A) annotation(
        Line(points = {{485, 145}, {510, 145}, {510, 143}}, color = {0, 0, 127}));
  connect(inversePark_fz.beta, clarkInv.B) annotation(
        Line(points = {{485, 137}, {510, 137}}, color = {0, 0, 127}));
  connect(constant2.y, clarkInv.C) annotation(
        Line(points = {{500, 98}, {503, 98}, {503, 132}, {510, 132}}, color = {0, 0, 127}));
      annotation(
        Diagram(coordinateSystem(extent = {{-150, 300}, {700, -130}}, grid = {1, 1})),
        Icon(coordinateSystem(extent = {{-120, -120}, {480, 160}}, grid = {1, 1})));
    end ART32;
    
    model ART52
      parameter Real Rf = 3.53e-4;
      parameter Real Lf = 2.81e-5;
      parameter Real Rr = 1.47;
      parameter Real Lr = 7.32e-3;
      parameter Real Rcc = 0.12;
      
      VSC_FZ.Testes.VSC vsc annotation(
        Placement(visible = true, transformation(origin = {-1, 27}, extent = {{-13, -13}, {13, 13}}, rotation = 0)));
      Modelica.Electrical.Analog.Sources.ConstantVoltage constantVoltage(V = 395.18) annotation(
        Placement(visible = true, transformation(origin = {-86, 18}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      Modelica.Electrical.Analog.Basic.Resistor resistor(R = Rf, alpha = 0) annotation(
        Placement(visible = true, transformation(origin = {87, 34}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Resistor resistor1(R = Rf, alpha = 0) annotation(
        Placement(visible = true, transformation(origin = {90, -26}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Resistor resistor2(R = Rf, alpha = 0) annotation(
        Placement(visible = true, transformation(origin = {92, -86}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Inductor inductor(L = Lf) annotation(
        Placement(visible = true, transformation(origin = {115, 34}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Inductor inductor1(L = Lf) annotation(
        Placement(visible = true, transformation(origin = {118, -26}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      VSC_FZ.Testes.Clark3p clark3p annotation(
        Placement(visible = true, transformation(origin = {-87, 129}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.VoltageSensor vsa annotation(
        Placement(visible = true, transformation(origin = {142, 53}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.VoltageSensor vsb annotation(
        Placement(visible = true, transformation(origin = {141, -6}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.VoltageSensor vscv annotation(
        Placement(visible = true, transformation(origin = {138, -64}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.CurrentSensor ia annotation(
        Placement(visible = true, transformation(origin = {54, 34}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.CurrentSensor ib annotation(
        Placement(visible = true, transformation(origin = {56, -26}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.CurrentSensor ic annotation(
        Placement(visible = true, transformation(origin = {58, -86}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
      VSC_FZ.Testes.Clark3p clark3p1 annotation(
        Placement(visible = true, transformation(origin = {205, 144}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression(y = ia.i) annotation(
        Placement(visible = true, transformation(origin = {171, 144}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression1(y = ib.i) annotation(
        Placement(visible = true, transformation(origin = {171, 128}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression2(y = ic.i) annotation(
        Placement(visible = true, transformation(origin = {171, 112}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression3(y = pLL_fz2.theta) annotation(
        Placement(visible = true, transformation(origin = {215.5, 116}, extent = {{-20.5, -9}, {20.5, 9}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression4(y = vsb.v) annotation(
        Placement(visible = true, transformation(origin = {-129, 123}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression5(y = vscv.v) annotation(
        Placement(visible = true, transformation(origin = {-129, 107}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression6(y = vsa.v) annotation(
        Placement(visible = true, transformation(origin = {-129, 139}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Sources.SineVoltage sineVoltage(V = 11.268e3, freqHz = 60) annotation(
        Placement(visible = true, transformation(origin = {362, 39}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Inductor inductor2(L = Lf) annotation(
        Placement(visible = true, transformation(origin = {120, -86}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Sources.SineVoltage sineVoltage1(V = 11.268e3, freqHz = 60, phase = -2.094395102393195) annotation(
        Placement(visible = true, transformation(origin = {364, -21}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Sources.SineVoltage sineVoltage2(V = 11.268e3, freqHz = 60, phase = 2.094395102393195) annotation(
        Placement(visible = true, transformation(origin = {364, -81}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Ground ground1 annotation(
        Placement(visible = true, transformation(origin = {399, -21}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
      Modelica.Blocks.Sources.RealExpression realExpression7(y = pLL_fz2.d) annotation(
        Placement(visible = true, transformation(origin = {344, 173.5}, extent = {{-17, -9.5}, {17, 9.5}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression8(y = pLL_fz2.q) annotation(
        Placement(visible = true, transformation(origin = {311.5, 98.5}, extent = {{-15.5, -9.5}, {15.5, 9.5}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression9(y = Controle_potencia.idref) annotation(
        Placement(visible = true, transformation(origin = {333, 158.5}, extent = {{-28, -9.5}, {28, 9.5}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression10(y = Correntes.d) annotation(
        Placement(visible = true, transformation(origin = {341, 145}, extent = {{-19, -10}, {19, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression11(y = Correntes.q) annotation(
        Placement(visible = true, transformation(origin = {341.5, 126}, extent = {{-19.5, -10}, {19.5, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression12(y = Controle_potencia.iqref) annotation(
        Placement(visible = true, transformation(origin = {335, 113.5}, extent = {{-26, -8.5}, {26, 8.5}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression13(y = pLL_fz2.theta) annotation(
        Placement(visible = true, transformation(origin = {445, 97.5}, extent = {{-23, -11.5}, {23, 11.5}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant2(k = 0) annotation(
        Placement(visible = true, transformation(origin = {489, 98}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression14(y = clarkInv.a) annotation(
        Placement(visible = true, transformation(origin = {-73, -56.5}, extent = {{-17, -8.5}, {17, 8.5}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression15(y = clarkInv.b) annotation(
        Placement(visible = true, transformation(origin = {-75, -83.5}, extent = {{-17, -9.5}, {17, 9.5}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression16(y = clarkInv.c) annotation(
        Placement(visible = true, transformation(origin = {-76, -111.5}, extent = {{-18, -9.5}, {18, 9.5}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.VoltageSensor vdc_2 annotation(
        Placement(visible = true, transformation(origin = {-109, 21}, extent = {{-7, -7}, {7, 7}}, rotation = -90)));
      Modelica.Blocks.Math.Division division annotation(
        Placement(visible = true, transformation(origin = {434, 157}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Math.Division division1 annotation(
        Placement(visible = true, transformation(origin = {432, 121}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Nonlinear.Limiter limiter(limitsAtInit = true, uMax = 1, uMin = -1) annotation(
        Placement(visible = true, transformation(origin = {-28, -58}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
      Modelica.Blocks.Nonlinear.Limiter limiter1(limitsAtInit = true, uMax = 1, uMin = -1) annotation(
        Placement(visible = true, transformation(origin = {-28, -84}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
      Modelica.Blocks.Nonlinear.Limiter limiter2(limitsAtInit = true, uMax = 1, uMin = -1) annotation(
        Placement(visible = true, transformation(origin = {-28, -112}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
      VSC_FZ.Testes.InnerControl2 Controle_corrente(Ictrl_KI = 0.5545, Ictrl_KP = 4.42e-2, lw = 2*3.1416*60*Lf, r = Rf) annotation(
        Placement(visible = true, transformation(origin = {394.5, 145.5}, extent = {{-14.5, -14.5}, {14.5, 14.5}}, rotation = 0)));
      Modelica.Blocks.Sources.Step Pref(height = 915e3, startTime = 0.1) annotation(
        Placement(visible = true, transformation(origin = {38, 139}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Step Qref(height = 0, startTime = 0) annotation(
        Placement(visible = true, transformation(origin = {37, 72}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      VSC_FZ.Testes.InversePark_fz inversePark_fz annotation(
        Placement(visible = true, transformation(origin = {474, 141}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      VSC_FZ.Testes.Park_Fz Correntes annotation(
        Placement(visible = true, transformation(origin = {243, 140}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      VSC_FZ.Testes.PLL_fz2 pLL_fz2(frequency = 60) annotation(
        Placement(visible = true, transformation(origin = {-34, 133}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Sources.ConstantVoltage constantVoltage1(V = 395.18) annotation(
        Placement(visible = true, transformation(origin = {-86, -16}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      Modelica.Electrical.Analog.Basic.Ground ground annotation(
        Placement(visible = true, transformation(origin = {-65, 70}, extent = {{-54, -86}, {-38, -70}}, rotation = 0)));
      Modelica.Blocks.Continuous.Filter filter(analogFilter = Modelica.Blocks.Types.AnalogFilter.Bessel, f_cut = 60, filterType = Modelica.Blocks.Types.FilterType.LowPass) annotation(
        Placement(visible = true, transformation(origin = {406, -18}, extent = {{84, 4}, {104, 24}}, rotation = 0)));
      VSC_FZ.Testes.PowerCalc powerCalc annotation(
        Placement(visible = true, transformation(origin = {465, -11}, extent = {{-14, -14}, {14, 14}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression17(y = vsa.v) annotation(
        Placement(visible = true, transformation(origin = {431, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression24(y = vsb.v) annotation(
        Placement(visible = true, transformation(origin = {431, -3}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression25(y = vscv.v) annotation(
        Placement(visible = true, transformation(origin = {431, -16}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression26(y = ia.i) annotation(
        Placement(visible = true, transformation(origin = {431, -31}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression27(y = ib.i) annotation(
        Placement(visible = true, transformation(origin = {431, -47}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression28(y = ic.i) annotation(
        Placement(visible = true, transformation(origin = {431, -63}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Continuous.Filter filter1(analogFilter = Modelica.Blocks.Types.AnalogFilter.Bessel, f_cut = 60, filterType = Modelica.Blocks.Types.FilterType.LowPass) annotation(
        Placement(visible = true, transformation(origin = {406, -49}, extent = {{84, 4}, {104, 24}}, rotation = 0)));
      VSC_FZ.Testes.PowerControl1 Controle_potencia(Pctrl_KI = 0.168, Pctrl_KP = 5.36e-4) annotation(
        Placement(visible = true, transformation(origin = {108.5, 105.5}, extent = {{-17.5, -17.5}, {17.5, 17.5}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression23(y = powerCalc.P) annotation(
        Placement(visible = true, transformation(origin = {33.5, 111}, extent = {{-20.5, -10}, {20.5, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression29(y = powerCalc.Q) annotation(
        Placement(visible = true, transformation(origin = {33.5, 94}, extent = {{-20.5, -10}, {20.5, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression18(y = vdc_2.v) annotation(
        Placement(visible = true, transformation(origin = {388, 78.5}, extent = {{-23, -11.5}, {23, 11.5}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Resistor resistor3(R = Rr, alpha = 0) annotation(
        Placement(visible = true, transformation(origin = {299, 39}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Resistor resistor4(R = Rr, alpha = 0) annotation(
        Placement(visible = true, transformation(origin = {298, -21}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Resistor resistor5(R = Rr, alpha = 0) annotation(
        Placement(visible = true, transformation(origin = {298, -81}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Ground ground2 annotation(
        Placement(visible = true, transformation(origin = {161.5, -6.5}, extent = {{-6.5, -6.5}, {6.5, 6.5}}, rotation = 90)));
      Modelica.Electrical.Analog.Basic.Ground ground3 annotation(
        Placement(visible = true, transformation(origin = {160.5, 52.5}, extent = {{-6.5, -6.5}, {6.5, 6.5}}, rotation = 90)));
      Modelica.Electrical.Analog.Basic.Ground ground4 annotation(
        Placement(visible = true, transformation(origin = {160.5, -64.5}, extent = {{-6.5, -6.5}, {6.5, 6.5}}, rotation = 90)));
      Modelica.Electrical.Analog.Basic.Inductor inductor3(L = Lr) annotation(
        Placement(visible = true, transformation(origin = {331, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Inductor inductor4(L = Lr) annotation(
        Placement(visible = true, transformation(origin = {330, -21}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Inductor inductor5(L = Lr) annotation(
        Placement(visible = true, transformation(origin = {330, -81}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      VSC_FZ.Testes.y_D_Transformer_ideal y_D_Transformer_ideal(Vac_primary = 440, Vac_secondary = 13800) annotation(
        Placement(visible = true, transformation(origin = {208, -28}, extent = {{-22, -22}, {22, 22}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Ground ground5 annotation(
        Placement(visible = true, transformation(origin = {302.5, 14.5}, extent = {{-6.5, -6.5}, {6.5, 6.5}}, rotation = 90)));
      Modelica.Electrical.Analog.Sensors.VoltageSensor pca annotation(
        Placement(visible = true, transformation(origin = {284, 15}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.VoltageSensor pcc annotation(
        Placement(visible = true, transformation(origin = {282, -106}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Ground ground6 annotation(
        Placement(visible = true, transformation(origin = {300.5, -106.5}, extent = {{-6.5, -6.5}, {6.5, 6.5}}, rotation = 90)));
      Modelica.Electrical.Analog.Basic.Ground ground7 annotation(
        Placement(visible = true, transformation(origin = {305.5, -44.5}, extent = {{-6.5, -6.5}, {6.5, 6.5}}, rotation = 90)));
      Modelica.Electrical.Analog.Sensors.VoltageSensor pcb annotation(
        Placement(visible = true, transformation(origin = {287, -44}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.CurrentSensor currentSensor annotation(
        Placement(visible = true, transformation(origin = {-97, 1}, extent = {{-4, -4}, {4, 4}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Ground ground8 annotation(
        Placement(visible = true, transformation(origin = {92.5, 9.5}, extent = {{-6.5, -6.5}, {6.5, 6.5}}, rotation = 90)));
      Modelica.Electrical.Analog.Sensors.VoltageSensor vta annotation(
        Placement(visible = true, transformation(origin = {74, 10}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.VoltageSensor vtb annotation(
        Placement(visible = true, transformation(origin = {78, -51}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Ground ground9 annotation(
        Placement(visible = true, transformation(origin = {96.5, -51.5}, extent = {{-6.5, -6.5}, {6.5, 6.5}}, rotation = 90)));
      Modelica.Electrical.Analog.Basic.Ground ground10 annotation(
        Placement(visible = true, transformation(origin = {101.5, -115.5}, extent = {{-6.5, -6.5}, {6.5, 6.5}}, rotation = 90)));
      Modelica.Electrical.Analog.Sensors.VoltageSensor vtc annotation(
        Placement(visible = true, transformation(origin = {83, -115}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Math.Division pcd annotation(
        Placement(visible = true, transformation(origin = {213, 210}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression30(y = pcc.v) annotation(
        Placement(visible = true, transformation(origin = {77, 190}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression31(y = pcb.v) annotation(
        Placement(visible = true, transformation(origin = {77, 206}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression32(y = pca.v) annotation(
        Placement(visible = true, transformation(origin = {77, 222}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Division pcq annotation(
        Placement(visible = true, transformation(origin = {212, 186}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant3(k = 11.267e3) annotation(
        Placement(visible = true, transformation(origin = {173, 185}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      VSC_FZ.Testes.Clark3p clark3p3 annotation(
        Placement(visible = true, transformation(origin = {121, 209}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      VSC_FZ.Testes.PLL_fz2 pLL_fz21(frequency = 60) annotation(
        Placement(visible = true, transformation(origin = {174, 213}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Division vsq annotation(
        Placement(visible = true, transformation(origin = {-16, 81}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Math.Division vsd annotation(
        Placement(visible = true, transformation(origin = {-15, 105}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant4(k = 359.258) annotation(
        Placement(visible = true, transformation(origin = {-55, 80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Division Pref_pu annotation(
        Placement(visible = true, transformation(origin = {130, 159}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant5(k = 915e3) annotation(
        Placement(visible = true, transformation(origin = {97, 143}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant6(k = 915e3) annotation(
        Placement(visible = true, transformation(origin = {506, 31}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Division P_pu annotation(
        Placement(visible = true, transformation(origin = {539, 47}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant7(k = 915e3) annotation(
        Placement(visible = true, transformation(origin = {506, -88}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Division Q_pu annotation(
        Placement(visible = true, transformation(origin = {539, -72}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant8(k = 1698) annotation(
        Placement(visible = true, transformation(origin = {263, 173}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Division id_pu annotation(
        Placement(visible = true, transformation(origin = {296, 189}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant9(k = 1698) annotation(
        Placement(visible = true, transformation(origin = {238, 95}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Division iq_pu annotation(
        Placement(visible = true, transformation(origin = {271, 111}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant10(k = 1698) annotation(
        Placement(visible = true, transformation(origin = {356, 205}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Division idref_pu annotation(
        Placement(visible = true, transformation(origin = {389, 221}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression19(y = Controle_potencia.idref) annotation(
        Placement(visible = true, transformation(origin = {326, 232.5}, extent = {{-28, -9.5}, {28, 9.5}}, rotation = 0)));
      Modelica.Blocks.Math.Division iqref_pu annotation(
        Placement(visible = true, transformation(origin = {503, 221}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant11(k = 1698) annotation(
        Placement(visible = true, transformation(origin = {470, 205}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression21(y = Controle_potencia.iqref) annotation(
        Placement(visible = true, transformation(origin = {444, 229.5}, extent = {{-26, -8.5}, {26, 8.5}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant12(k = 915e3) annotation(
        Placement(visible = true, transformation(origin = {116, 259}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Division Qref_pu annotation(
        Placement(visible = true, transformation(origin = {149, 275}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression20(y = Qref.y) annotation(
        Placement(visible = true, transformation(origin = {107, 283.5}, extent = {{-17, -9.5}, {17, 9.5}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.CurrentSensor ipca annotation(
        Placement(visible = true, transformation(origin = {260, 39}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.CurrentSensor ipcb annotation(
        Placement(visible = true, transformation(origin = {262, -21}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.CurrentSensor ipcc annotation(
        Placement(visible = true, transformation(origin = {264, -81}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression39(y = ia.i) annotation(
        Placement(visible = true, transformation(origin = {-118, 245}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression40(y = ib.i) annotation(
        Placement(visible = true, transformation(origin = {-118, 229}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      VSC_FZ.Testes.PLL_fz2 correntes_tot(frequency = 60) annotation(
        Placement(visible = true, transformation(origin = {-18, 243}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant14(k = 1698) annotation(
        Placement(visible = true, transformation(origin = {-39, 190}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Division itot_d annotation(
        Placement(visible = true, transformation(origin = {1, 215}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Math.Division itot_q annotation(
        Placement(visible = true, transformation(origin = {0, 191}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      VSC_FZ.Testes.Clark3p clark3p2 annotation(
        Placement(visible = true, transformation(origin = {-71, 239}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression38(y = ic.i) annotation(
        Placement(visible = true, transformation(origin = {-118, 213}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant1(k = 915e3) annotation(
        Placement(visible = true, transformation(origin = {655, 34}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression22(y = ipcb.i) annotation(
        Placement(visible = true, transformation(origin = {580, -44}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant13(k = 915e3) annotation(
        Placement(visible = true, transformation(origin = {655, -85}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Continuous.Filter filter2(analogFilter = Modelica.Blocks.Types.AnalogFilter.Bessel, f_cut = 60, filterType = Modelica.Blocks.Types.FilterType.LowPass) annotation(
        Placement(visible = true, transformation(origin = {555, -46}, extent = {{84, 4}, {104, 24}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression33(y = ipca.i) annotation(
        Placement(visible = true, transformation(origin = {580, -28}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Division P_pupc annotation(
        Placement(visible = true, transformation(origin = {688, 50}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Continuous.Filter filter3(analogFilter = Modelica.Blocks.Types.AnalogFilter.Bessel, f_cut = 60, filterType = Modelica.Blocks.Types.FilterType.LowPass) annotation(
        Placement(visible = true, transformation(origin = {555, -15}, extent = {{84, 4}, {104, 24}}, rotation = 0)));
      Modelica.Blocks.Math.Division Q_pupc annotation(
        Placement(visible = true, transformation(origin = {688, -69}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression34(y = pcb.v) annotation(
        Placement(visible = true, transformation(origin = {580, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression35(y = ipcc.i) annotation(
        Placement(visible = true, transformation(origin = {580, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression36(y = pcc.v) annotation(
        Placement(visible = true, transformation(origin = {580, -13}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      VSC_FZ.Testes.PowerCalc powerCalc_pcc annotation(
        Placement(visible = true, transformation(origin = {614, -8}, extent = {{-14, -14}, {14, 14}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression37(y = pca.v) annotation(
        Placement(visible = true, transformation(origin = {580, 13}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      VSC_FZ.Testes.PLL_fz2 ipc(frequency = 60) annotation(
        Placement(visible = true, transformation(origin = {611, 276}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression41(y = ipcc.i) annotation(
        Placement(visible = true, transformation(origin = {516, 243}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Division ipd annotation(
        Placement(visible = true, transformation(origin = {677, 243}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression42(y = ipcb.i) annotation(
        Placement(visible = true, transformation(origin = {511, 262}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      VSC_FZ.Testes.Clark3p clark3p4 annotation(
        Placement(visible = true, transformation(origin = {560, 274}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant15(k = 54.13) annotation(
        Placement(visible = true, transformation(origin = {630, 217}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression43(y = ipca.i) annotation(
        Placement(visible = true, transformation(origin = {511, 279}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Division ipq annotation(
        Placement(visible = true, transformation(origin = {669, 218}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Electrical.Analog.Ideal.IdealClosingSwitch idealClosingSwitch annotation(
        Placement(visible = true, transformation(origin = {337, -47}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
      Modelica.Electrical.Analog.Basic.Resistor resistor6(R = Rcc, alpha = 0) annotation(
        Placement(visible = true, transformation(origin = {352, 1}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Ground ground11 annotation(
        Placement(visible = true, transformation(origin = {372.5, 1.5}, extent = {{-6.5, -6.5}, {6.5, 6.5}}, rotation = 90)));
      Modelica.Electrical.Analog.Ideal.IdealClosingSwitch idealClosingSwitch1 annotation(
        Placement(visible = true, transformation(origin = {337, -105}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
      Modelica.Electrical.Analog.Ideal.IdealClosingSwitch switch2 annotation(
        Placement(visible = true, transformation(origin = {334, 11}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
      Modelica.Blocks.MathBoolean.And and1(nu = 2) annotation(
        Placement(visible = true, transformation(origin = {294, -155}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
      Modelica.Blocks.Sources.BooleanStep booleanStep4(startTime = 0.2082, startValue = false) annotation(
        Placement(visible = true, transformation(origin = {268, -140}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Ground ground13 annotation(
        Placement(visible = true, transformation(origin = {374.5, -114.5}, extent = {{-6.5, -6.5}, {6.5, 6.5}}, rotation = 90)));
      Modelica.Electrical.Analog.Basic.Resistor resistor7(R = Rcc, alpha = 0) annotation(
        Placement(visible = true, transformation(origin = {354, -57}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Ground ground12 annotation(
        Placement(visible = true, transformation(origin = {374.5, -56.5}, extent = {{-6.5, -6.5}, {6.5, 6.5}}, rotation = 90)));
      Modelica.Blocks.Sources.BooleanStep booleanStep3(startTime = 0, startValue = false) annotation(
        Placement(visible = true, transformation(origin = {266, -173}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Resistor resistor8(R = Rcc, alpha = 0) annotation(
        Placement(visible = true, transformation(origin = {354, -117}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.BooleanStep booleanStep(startTime = 0.21385, startValue = false) annotation(
        Placement(visible = true, transformation(origin = {194, -140}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.MathBoolean.And and2(nu = 2) annotation(
        Placement(visible = true, transformation(origin = {220, -155}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
    Modelica.Blocks.Sources.BooleanStep booleanStep1(startTime = 0, startValue = false) annotation(
        Placement(visible = true, transformation(origin = {192, -173}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.BooleanStep booleanStep2(startTime = 0.21108, startValue = false) annotation(
        Placement(visible = true, transformation(origin = {346, -157}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.MathBoolean.And and3(nu = 2) annotation(
        Placement(visible = true, transformation(origin = {372, -172}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
    Modelica.Blocks.Sources.BooleanStep booleanStep5(startTime = 0, startValue = false) annotation(
        Placement(visible = true, transformation(origin = {344, -190}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Division ib_pu annotation(
        Placement(visible = true, transformation(origin = {276, 286}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
  Modelica.Blocks.Sources.RealExpression realExpression46(y = ia.i) annotation(
        Placement(visible = true, transformation(origin = {219, 304}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.RealExpression realExpression44(y = ib.i) annotation(
        Placement(visible = true, transformation(origin = {219, 288}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.RealExpression realExpression45(y = ic.i) annotation(
        Placement(visible = true, transformation(origin = {219, 272}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Constant constant16(k = 1698) annotation(
        Placement(visible = true, transformation(origin = {239, 241}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Division ia_pu annotation(
        Placement(visible = true, transformation(origin = {276, 313}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
  Modelica.Blocks.Math.Division ic_pu annotation(
        Placement(visible = true, transformation(origin = {276, 260}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
  ClarkInv clarkInv annotation(
        Placement(visible = true, transformation(origin = {530, 141}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(constantVoltage.p, vsc.pin_p) annotation(
        Line(points = {{-86, 28}, {-86, 34}, {-15, 34}, {-15, 33.5}}, color = {0, 0, 255}));
      connect(resistor.n, inductor.p) annotation(
        Line(points = {{97, 34}, {105, 34}}, color = {0, 0, 255}));
      connect(resistor1.n, inductor1.p) annotation(
        Line(points = {{100, -26}, {108, -26}}, color = {0, 0, 255}));
      connect(vsc.Vta, ia.p) annotation(
        Line(points = {{13, 35}, {30, 35}, {30, 34}, {48, 34}}, color = {0, 0, 255}));
      connect(ia.n, resistor.p) annotation(
        Line(points = {{60, 34}, {77, 34}}, color = {0, 0, 255}));
      connect(vsc.Vtb, ib.p) annotation(
        Line(points = {{13, 27}, {40, 27}, {40, -26}, {50, -26}}, color = {0, 0, 255}));
      connect(ib.n, resistor1.p) annotation(
        Line(points = {{62, -26}, {80, -26}}, color = {0, 0, 255}));
      connect(ic.n, resistor2.p) annotation(
        Line(points = {{64, -86}, {82, -86}}, color = {0, 0, 255}));
      connect(ic.p, vsc.Vtc) annotation(
        Line(points = {{52, -86}, {24, -86}, {24, 19}, {13, 19}}, color = {0, 0, 255}));
      connect(realExpression.y, clark3p1.V_abc_D[1]) annotation(
        Line(points = {{182, 144}, {192, 144}}, color = {0, 0, 127}));
      connect(realExpression1.y, clark3p1.V_abc_D[2]) annotation(
        Line(points = {{182, 128}, {184, 128}, {184, 144}, {192, 144}}, color = {0, 0, 127}));
      connect(realExpression2.y, clark3p1.V_abc_D[3]) annotation(
        Line(points = {{182, 112}, {184, 112}, {184, 144}, {192, 144}}, color = {0, 0, 127}));
      connect(realExpression6.y, clark3p.V_abc_D[1]) annotation(
        Line(points = {{-118, 139}, {-104, 139}, {-104, 129}, {-99, 129}}, color = {0, 0, 127}));
      connect(realExpression4.y, clark3p.V_abc_D[2]) annotation(
        Line(points = {{-118, 123}, {-109, 123}, {-109, 129}, {-99, 129}}, color = {0, 0, 127}));
      connect(realExpression5.y, clark3p.V_abc_D[3]) annotation(
        Line(points = {{-118, 107}, {-112, 107}, {-112, 129}, {-99, 129}}, color = {0, 0, 127}));
      connect(inductor2.p, resistor2.n) annotation(
        Line(points = {{110, -86}, {102, -86}}, color = {0, 0, 255}));
      connect(sineVoltage.n, ground1.p) annotation(
        Line(points = {{372, 39}, {389, 39}, {389, -21}}, color = {0, 0, 255}));
      connect(sineVoltage1.n, ground1.p) annotation(
        Line(points = {{374, -21}, {389, -21}}, color = {0, 0, 255}));
      connect(sineVoltage2.n, ground1.p) annotation(
        Line(points = {{374, -81}, {389, -81}, {389, -21}}, color = {0, 0, 255}));
      connect(vsa.p, inductor.n) annotation(
        Line(points = {{135, 53}, {125, 53}, {125, 34}}, color = {0, 0, 255}));
      connect(vsb.p, inductor1.n) annotation(
        Line(points = {{134, -6}, {128, -6}, {128, -26}}, color = {0, 0, 255}));
      connect(vscv.p, inductor2.n) annotation(
        Line(points = {{130, -64}, {130, -86}}, color = {0, 0, 255}));
      connect(vdc_2.p, constantVoltage.p) annotation(
        Line(points = {{-109, 28}, {-87, 28}}, color = {0, 0, 255}));
      connect(vdc_2.n, constantVoltage.n) annotation(
        Line(points = {{-109, 14}, {-109, 8}, {-87, 8}}, color = {0, 0, 255}));
      connect(realExpression14.y, limiter.u) annotation(
        Line(points = {{-54, -56.5}, {-46, -56.5}, {-46, -58}, {-38, -58}}, color = {0, 0, 127}));
      connect(limiter.y, vsc.ma) annotation(
        Line(points = {{-19, -58}, {-8, -58}, {-8, 12}}, color = {0, 0, 127}));
      connect(realExpression15.y, limiter1.u) annotation(
        Line(points = {{-56, -83.5}, {-47, -83.5}, {-47, -84}, {-38, -84}}, color = {0, 0, 127}));
      connect(realExpression16.y, limiter2.u) annotation(
        Line(points = {{-56, -111.5}, {-47, -111.5}, {-47, -112}, {-38, -112}}, color = {0, 0, 127}));
      connect(limiter1.y, vsc.mb) annotation(
        Line(points = {{-20, -84}, {0, -84}, {0, 12}}, color = {0, 0, 127}));
      connect(limiter2.y, vsc.mc) annotation(
        Line(points = {{-20, -112}, {6, -112}, {6, 12}}, color = {0, 0, 127}));
      connect(Controle_corrente.vdref, division.u1) annotation(
        Line(points = {{410, 154}, {410, 161.22}, {424.45, 161.22}}, color = {0, 0, 127}));
      connect(Controle_corrente.vqref, division1.u1) annotation(
        Line(points = {{410, 139}, {410, 125.23}, {422.45, 125.23}}, color = {0, 0, 127}));
      connect(realExpression7.y, Controle_corrente.vd) annotation(
        Line(points = {{362.7, 173.5}, {373.7, 173.5}, {373.7, 159}, {379, 159}}, color = {0, 0, 127}));
      connect(realExpression10.y, Controle_corrente.Id) annotation(
        Line(points = {{361.9, 145}, {370.4, 145}, {370.4, 147}, {379, 147}}, color = {0, 0, 127}));
      connect(realExpression11.y, Controle_corrente.Iq) annotation(
        Line(points = {{362.95, 126}, {369.95, 126}, {369.95, 142}, {379, 142}}, color = {0, 0, 127}));
      connect(realExpression8.y, Controle_corrente.vq) annotation(
        Line(points = {{328.55, 98.5}, {328.55, 99.5}, {379, 99.5}, {379, 132}}, color = {0, 0, 127}));
      connect(division.y, inversePark_fz.d) annotation(
        Line(points = {{441.7, 157}, {453.7, 157}, {453.7, 145}, {461.7, 145}}, color = {0, 0, 127}));
      connect(division1.y, inversePark_fz.q) annotation(
        Line(points = {{439.7, 121}, {455.7, 121}, {455.7, 137}, {461.7, 137}}, color = {0, 0, 127}));
      connect(inversePark_fz.theta, realExpression13.y) annotation(
        Line(points = {{474, 129}, {474, 97.5}, {470, 97.5}}, color = {0, 0, 127}));
      connect(clark3p1.V_ab_D[1], Correntes.alpha) annotation(
        Line(points = {{216, 144}, {231, 144}}, color = {0, 0, 127}));
      connect(clark3p1.V_ab_D[2], Correntes.beta) annotation(
        Line(points = {{216, 144}, {224, 144}, {224, 136}, {231, 136}}, color = {0, 0, 127}));
      connect(Correntes.theta, realExpression3.y) annotation(
        Line(points = {{243, 128}, {243, 116}, {238, 116}}, color = {0, 0, 127}));
      connect(clark3p.V_ab_D[1], pLL_fz2.alpha) annotation(
        Line(points = {{-76, 129}, {-63.5, 129}, {-63.5, 137}, {-46, 137}}, color = {0, 0, 127}));
      connect(pLL_fz2.beta, clark3p.V_ab_D[2]) annotation(
        Line(points = {{-46.2, 128.4}, {-59.2, 128.4}, {-59.2, 129}, {-76, 129}}, color = {0, 0, 127}));
      connect(constantVoltage1.n, vsc.pin_n) annotation(
        Line(points = {{-86, -26}, {-60, -26}, {-60, 20.5}, {-15, 20.5}}, color = {0, 0, 255}));
      connect(constantVoltage.n, constantVoltage1.p) annotation(
        Line(points = {{-86, 8}, {-86, -6}}, color = {0, 0, 255}));
      connect(realExpression17.y, powerCalc.V1) annotation(
        Line(points = {{442, 10}, {447, 10}, {447, 1}, {450, 1}}, color = {0, 0, 127}));
      connect(realExpression24.y, powerCalc.V2) annotation(
        Line(points = {{442, -3}, {450, -3}}, color = {0, 0, 127}));
      connect(realExpression25.y, powerCalc.V3) annotation(
        Line(points = {{442, -16}, {444, -16}, {444, -7}, {450, -7}}, color = {0, 0, 127}));
      connect(realExpression26.y, powerCalc.I1) annotation(
        Line(points = {{442, -31}, {446, -31}, {446, -14}, {450, -14}}, color = {0, 0, 127}));
      connect(realExpression27.y, powerCalc.I2) annotation(
        Line(points = {{442, -47}, {447, -47}, {447, -18}, {450, -18}}, color = {0, 0, 127}));
      connect(realExpression28.y, powerCalc.I3) annotation(
        Line(points = {{442, -63}, {448, -63}, {448, -23}, {450, -23}}, color = {0, 0, 127}));
      connect(powerCalc.P, filter.u) annotation(
        Line(points = {{480.4, -4.56}, {488.4, -4.56}, {488.4, -3.56}}, color = {0, 0, 127}));
      connect(powerCalc.Q, filter1.u) annotation(
        Line(points = {{480.4, -15.2}, {483.4, -15.2}, {483.4, -35.2}, {488.4, -35.2}}, color = {0, 0, 127}));
      connect(Qref.y, Controle_potencia.Q_ref) annotation(
        Line(points = {{48, 72}, {78, 72}, {78, 95}, {89, 95}}, color = {0, 0, 127}));
      connect(Pref.y, Controle_potencia.P_ref) annotation(
        Line(points = {{49, 139}, {81, 139}, {81, 114}, {89, 114}}, color = {0, 0, 127}));
      connect(realExpression9.y, Controle_corrente.Id_ref) annotation(
        Line(points = {{364, 158.5}, {369.25, 158.5}, {369.25, 153}, {379, 153}}, color = {0, 0, 127}));
      connect(realExpression12.y, Controle_corrente.Iqref) annotation(
        Line(points = {{364, 113.5}, {374.1, 113.5}, {374.1, 137}, {379, 137}}, color = {0, 0, 127}));
      connect(realExpression29.y, Controle_potencia.Q) annotation(
        Line(points = {{56, 94}, {72, 94}, {72, 101}, {89, 101}}, color = {0, 0, 127}));
      connect(realExpression23.y, Controle_potencia.P) annotation(
        Line(points = {{56, 111}, {67, 111}, {67, 107}, {89, 107}}, color = {0, 0, 127}));
      connect(realExpression18.y, division1.u2) annotation(
        Line(points = {{413.3, 78.5}, {413.3, 117}, {424.3, 117}}, color = {0, 0, 127}));
      connect(realExpression18.y, division.u2) annotation(
        Line(points = {{413.3, 78.5}, {413.3, 153}, {426.3, 153}}, color = {0, 0, 127}));
      connect(vsb.n, ground2.p) annotation(
        Line(points = {{148, -6}, {155, -6}}, color = {0, 0, 255}));
      connect(vsa.n, ground3.p) annotation(
        Line(points = {{149, 53}, {154, 53}}, color = {0, 0, 255}));
      connect(vscv.n, ground4.p) annotation(
        Line(points = {{146, -64}, {154, -64}}, color = {0, 0, 255}));
      connect(resistor4.n, inductor4.n) annotation(
        Line(points = {{308, -21}, {320, -21}}, color = {0, 0, 255}));
      connect(resistor5.n, inductor5.n) annotation(
        Line(points = {{308, -81}, {320, -81}}, color = {0, 0, 255}));
      connect(inductor5.p, sineVoltage2.p) annotation(
        Line(points = {{340, -81}, {354, -81}}, color = {0, 0, 255}));
      connect(inductor4.p, sineVoltage1.p) annotation(
        Line(points = {{340, -21}, {354, -21}}, color = {0, 0, 255}));
      connect(resistor3.n, inductor3.p) annotation(
        Line(points = {{309, 39}, {313, 39}, {313, 40}, {321, 40}}, color = {0, 0, 255}));
      connect(inductor3.n, sineVoltage.p) annotation(
        Line(points = {{341, 40}, {352, 40}, {352, 39}}, color = {0, 0, 255}));
      connect(y_D_Transformer_ideal.pin_p, inductor.n) annotation(
        Line(points = {{184, -14}, {180, -14}, {180, 34}, {125, 34}}, color = {0, 0, 255}));
      connect(y_D_Transformer_ideal.pin_p1, inductor1.n) annotation(
        Line(points = {{184, -25}, {128, -25}, {128, -26}}, color = {0, 0, 255}));
      connect(y_D_Transformer_ideal.pin_p2, inductor2.n) annotation(
        Line(points = {{184, -35}, {171, -35}, {171, -86}, {130, -86}}, color = {0, 0, 255}));
      connect(pca.n, ground5.p) annotation(
        Line(points = {{291, 15}, {296, 15}}, color = {0, 0, 255}));
      connect(pcc.n, ground6.p) annotation(
        Line(points = {{289, -106}, {294, -106}}, color = {0, 0, 255}));
      connect(pcb.n, ground7.p) annotation(
        Line(points = {{294, -44}, {299, -44}}, color = {0, 0, 255}));
      connect(pca.p, y_D_Transformer_ideal.pin_n) annotation(
        Line(points = {{277, 15}, {238, 15}, {238, -15}, {232, -15}}, color = {0, 0, 255}));
      connect(pcb.p, y_D_Transformer_ideal.pin_n1) annotation(
        Line(points = {{280, -44}, {248, -44}, {248, -24}, {232, -24}}, color = {0, 0, 255}));
      connect(pcc.p, y_D_Transformer_ideal.pin_n2) annotation(
        Line(points = {{275, -106}, {243, -106}, {243, -34}, {232, -34}}, color = {0, 0, 255}));
      connect(ground.p, currentSensor.p) annotation(
        Line(points = {{-111, 0}, {-101, 0}, {-101, 1}}, color = {0, 0, 255}));
      connect(currentSensor.n, constantVoltage.n) annotation(
        Line(points = {{-93, 1}, {-86, 1}, {-86, 8}}, color = {0, 0, 255}));
      connect(vta.n, ground8.p) annotation(
        Line(points = {{81, 10}, {86, 10}}, color = {0, 0, 255}));
      connect(vtb.n, ground9.p) annotation(
        Line(points = {{85, -51}, {90, -51}}, color = {0, 0, 255}));
      connect(vtc.n, ground10.p) annotation(
        Line(points = {{90, -115}, {95, -115}}, color = {0, 0, 255}));
      connect(vtc.p, resistor2.p) annotation(
        Line(points = {{76, -115}, {70, -115}, {70, -86}, {82, -86}}, color = {0, 0, 255}));
      connect(vtb.p, resistor1.p) annotation(
        Line(points = {{71, -51}, {65, -51}, {65, -26}, {80, -26}}, color = {0, 0, 255}));
      connect(vta.p, resistor.p) annotation(
        Line(points = {{67, 10}, {64, 10}, {64, 34}, {77, 34}}, color = {0, 0, 255}));
      connect(realExpression31.y, clark3p3.V_abc_D[2]) annotation(
        Line(points = {{88, 206}, {94, 206}, {94, 209}, {109, 209}}, color = {0, 0, 127}));
      connect(pLL_fz21.q, pcq.u1) annotation(
        Line(points = {{185, 206.6}, {193, 206.6}, {193, 189.6}, {204, 189.6}}, color = {0, 0, 127}));
      connect(realExpression32.y, clark3p3.V_abc_D[1]) annotation(
        Line(points = {{88, 222}, {102, 222}, {102, 209}, {109, 209}}, color = {0, 0, 127}));
      connect(realExpression30.y, clark3p3.V_abc_D[3]) annotation(
        Line(points = {{88, 190}, {94, 190}, {94, 209}, {109, 209}}, color = {0, 0, 127}));
      connect(pLL_fz21.d, pcd.u1) annotation(
        Line(points = {{185, 213}, {205, 213}, {205, 214}}, color = {0, 0, 127}));
      connect(clark3p3.V_ab_D[1], pLL_fz21.alpha) annotation(
        Line(points = {{132, 209}, {144.5, 209}, {144.5, 217}, {162, 217}}, color = {0, 0, 127}));
      connect(constant3.y, pcq.u2) annotation(
        Line(points = {{184, 185}, {204, 185}, {204, 182}}, color = {0, 0, 127}));
      connect(pLL_fz21.beta, clark3p3.V_ab_D[2]) annotation(
        Line(points = {{161.8, 208.4}, {148.8, 208.4}, {148.8, 209.4}, {131.8, 209.4}}, color = {0, 0, 127}));
      connect(constant3.y, pcd.u2) annotation(
        Line(points = {{184, 185}, {197, 185}, {197, 206}, {205, 206}}, color = {0, 0, 127}));
      connect(constant4.y, vsq.u2) annotation(
        Line(points = {{-44, 80}, {-24, 80}, {-24, 77}}, color = {0, 0, 127}));
      connect(constant4.y, vsd.u2) annotation(
        Line(points = {{-44, 80}, {-31, 80}, {-31, 101}, {-23, 101}}, color = {0, 0, 127}));
      connect(pLL_fz2.d, vsd.u1) annotation(
        Line(points = {{-23, 133}, {-33, 133}, {-33, 109}, {-23, 109}}, color = {0, 0, 127}));
      connect(pLL_fz2.q, vsq.u1) annotation(
        Line(points = {{-23, 127}, {-37, 127}, {-37, 85}, {-24, 85}}, color = {0, 0, 127}));
      connect(Pref_pu.u1, Pref.y) annotation(
        Line(points = {{122, 163}, {81, 163}, {81, 139}, {49, 139}}, color = {0, 0, 127}));
      connect(constant5.y, Pref_pu.u2) annotation(
        Line(points = {{108, 143}, {113, 143}, {113, 155}, {122, 155}}, color = {0, 0, 127}));
      connect(constant6.y, P_pu.u2) annotation(
        Line(points = {{517, 31}, {522, 31}, {522, 43}, {531, 43}}, color = {0, 0, 127}));
      connect(powerCalc.P, P_pu.u1) annotation(
        Line(points = {{480, -5}, {482, -5}, {482, 51}, {531, 51}}, color = {0, 0, 127}));
      connect(constant7.y, Q_pu.u2) annotation(
        Line(points = {{517, -88}, {522, -88}, {522, -76}, {531, -76}}, color = {0, 0, 127}));
      connect(powerCalc.Q, Q_pu.u1) annotation(
        Line(points = {{480, -15}, {484, -15}, {484, -68}, {531, -68}}, color = {0, 0, 127}));
      connect(constant8.y, id_pu.u2) annotation(
        Line(points = {{274, 173}, {279, 173}, {279, 185}, {288, 185}}, color = {0, 0, 127}));
      connect(id_pu.u1, Correntes.d) annotation(
        Line(points = {{288, 193}, {242, 193}, {242, 152}, {263, 152}, {263, 144}, {254, 144}}, color = {0, 0, 127}));
      connect(constant9.y, iq_pu.u2) annotation(
        Line(points = {{249, 95}, {254, 95}, {254, 107}, {263, 107}}, color = {0, 0, 127}));
      connect(iq_pu.u1, Correntes.q) annotation(
        Line(points = {{263, 115}, {257, 115}, {257, 136}, {254, 136}}, color = {0, 0, 127}));
      connect(constant10.y, idref_pu.u2) annotation(
        Line(points = {{367, 205}, {372, 205}, {372, 217}, {381, 217}}, color = {0, 0, 127}));
      connect(realExpression19.y, idref_pu.u1) annotation(
        Line(points = {{357, 233}, {375, 233}, {375, 225}, {381, 225}}, color = {0, 0, 127}));
      connect(constant11.y, iqref_pu.u2) annotation(
        Line(points = {{481, 205}, {486, 205}, {486, 217}, {495, 217}}, color = {0, 0, 127}));
      connect(realExpression21.y, iqref_pu.u1) annotation(
        Line(points = {{473, 229.5}, {484, 229.5}, {484, 225}, {495, 225}}, color = {0, 0, 127}));
      connect(constant12.y, Qref_pu.u2) annotation(
        Line(points = {{127, 259}, {132, 259}, {132, 271}, {141, 271}}, color = {0, 0, 127}));
      connect(realExpression20.y, Qref_pu.u1) annotation(
        Line(points = {{126, 284}, {133, 284}, {133, 279}, {141, 279}}, color = {0, 0, 127}));
      connect(resistor3.p, ipca.n) annotation(
        Line(points = {{289, 39}, {266, 39}}, color = {0, 0, 255}));
      connect(ipca.p, y_D_Transformer_ideal.pin_n) annotation(
        Line(points = {{254, 39}, {238, 39}, {238, -15}, {232, -15}}, color = {0, 0, 255}));
      connect(resistor4.p, ipcb.n) annotation(
        Line(points = {{288, -21}, {268, -21}}, color = {0, 0, 255}));
      connect(ipcb.p, y_D_Transformer_ideal.pin_n1) annotation(
        Line(points = {{256, -21}, {248, -21}, {248, -24}, {232, -24}}, color = {0, 0, 255}));
      connect(resistor5.p, ipcc.n) annotation(
        Line(points = {{288, -81}, {270, -81}}, color = {0, 0, 255}));
      connect(ipcc.p, y_D_Transformer_ideal.pin_n2) annotation(
        Line(points = {{258, -81}, {243, -81}, {243, -34}, {232, -34}}, color = {0, 0, 255}));
      connect(constant14.y, itot_q.u2) annotation(
        Line(points = {{-28, 190}, {-8, 190}, {-8, 187}}, color = {0, 0, 127}));
      connect(realExpression39.y, clark3p2.V_abc_D[1]) annotation(
        Line(points = {{-107, 245}, {-91, 245}, {-91, 239}, {-83, 239}}, color = {0, 0, 127}));
      connect(clark3p2.V_ab_D[1], correntes_tot.alpha) annotation(
        Line(points = {{-60, 239}, {-47.5, 239}, {-47.5, 247}, {-30, 247}}, color = {0, 0, 127}));
      connect(correntes_tot.q, itot_q.u1) annotation(
        Line(points = {{-7, 236.6}, {-21, 236.6}, {-21, 195}, {-8, 195}}, color = {0, 0, 127}));
      connect(constant14.y, itot_d.u2) annotation(
        Line(points = {{-28, 190}, {-15, 190}, {-15, 211}, {-7, 211}}, color = {0, 0, 127}));
      connect(realExpression40.y, clark3p2.V_abc_D[2]) annotation(
        Line(points = {{-107, 229}, {-91, 229}, {-91, 239}, {-83, 239}}, color = {0, 0, 127}));
      connect(correntes_tot.beta, clark3p2.V_ab_D[2]) annotation(
        Line(points = {{-30.2, 238.4}, {-43.2, 238.4}, {-43.2, 239.4}, {-60.2, 239.4}}, color = {0, 0, 127}));
      connect(realExpression38.y, clark3p2.V_abc_D[3]) annotation(
        Line(points = {{-107, 213}, {-91, 213}, {-91, 239}, {-83, 239}}, color = {0, 0, 127}));
      connect(correntes_tot.d, itot_d.u1) annotation(
        Line(points = {{-7, 243}, {-17, 243}, {-17, 219}, {-7, 219}}, color = {0, 0, 127}));
      connect(constant1.y, P_pupc.u2) annotation(
        Line(points = {{666, 34}, {671, 34}, {671, 46}, {680, 46}}, color = {0, 0, 127}));
      connect(realExpression33.y, powerCalc_pcc.I1) annotation(
        Line(points = {{591, -28}, {595, -28}, {595, -11}, {599, -11}}, color = {0, 0, 127}));
      connect(realExpression35.y, powerCalc_pcc.I3) annotation(
        Line(points = {{591, -60}, {597, -60}, {597, -20}, {599, -20}}, color = {0, 0, 127}));
      connect(constant13.y, Q_pupc.u2) annotation(
        Line(points = {{666, -85}, {671, -85}, {671, -73}, {680, -73}}, color = {0, 0, 127}));
      connect(realExpression37.y, powerCalc_pcc.V1) annotation(
        Line(points = {{591, 13}, {596, 13}, {596, 4}, {599, 4}}, color = {0, 0, 127}));
      connect(powerCalc_pcc.Q, Q_pupc.u1) annotation(
        Line(points = {{629.4, -12.2}, {633.4, -12.2}, {633.4, -65.2}, {680.4, -65.2}}, color = {0, 0, 127}));
      connect(realExpression22.y, powerCalc_pcc.I2) annotation(
        Line(points = {{591, -44}, {596, -44}, {596, -15}, {599, -15}}, color = {0, 0, 127}));
      connect(realExpression34.y, powerCalc_pcc.V2) annotation(
        Line(points = {{591, 0}, {599, 0}}, color = {0, 0, 127}));
      connect(realExpression36.y, powerCalc_pcc.V3) annotation(
        Line(points = {{591, -13}, {593, -13}, {593, -4}, {599, -4}}, color = {0, 0, 127}));
      connect(powerCalc_pcc.Q, filter2.u) annotation(
        Line(points = {{629.4, -12.2}, {632.4, -12.2}, {632.4, -32.2}, {637.4, -32.2}}, color = {0, 0, 127}));
      connect(powerCalc_pcc.P, filter3.u) annotation(
        Line(points = {{629.4, -1.56}, {637.4, -1.56}, {637.4, -0.56}}, color = {0, 0, 127}));
      connect(powerCalc_pcc.P, P_pupc.u1) annotation(
        Line(points = {{629.4, -1.56}, {631.4, -1.56}, {631.4, 54.44}, {680.4, 54.44}}, color = {0, 0, 127}));
      connect(realExpression43.y, clark3p4.V_abc_D[1]) annotation(
        Line(points = {{522, 279}, {541, 279}, {541, 274}, {548, 274}}, color = {0, 0, 127}));
      connect(ipc.q, ipq.u1) annotation(
        Line(points = {{622, 270}, {648, 270}, {648, 221.6}, {661, 221.6}}, color = {0, 0, 127}));
      connect(constant15.y, ipq.u2) annotation(
        Line(points = {{641, 217}, {661, 217}, {661, 214}}, color = {0, 0, 127}));
      connect(realExpression42.y, clark3p4.V_abc_D[2]) annotation(
        Line(points = {{522, 262}, {538, 262}, {538, 274}, {548, 274}}, color = {0, 0, 127}));
      connect(clark3p4.V_ab_D[1], ipc.alpha) annotation(
        Line(points = {{571, 274}, {581.5, 274}, {581.5, 280}, {599, 280}}, color = {0, 0, 127}));
      connect(realExpression41.y, clark3p4.V_abc_D[3]) annotation(
        Line(points = {{527, 243}, {539, 243}, {539, 274}, {548, 274}}, color = {0, 0, 127}));
      connect(ipc.d, ipd.u1) annotation(
        Line(points = {{622, 276}, {652, 276}, {652, 248}, {667, 248}}, color = {0, 0, 127}));
      connect(constant15.y, ipd.u2) annotation(
        Line(points = {{641, 217}, {654, 217}, {654, 238}, {667, 238}}, color = {0, 0, 127}));
      connect(ipc.beta, clark3p4.V_ab_D[2]) annotation(
        Line(points = {{598.8, 271.4}, {585.6, 271.4}, {585.6, 274}, {571, 274}}, color = {0, 0, 127}));
      connect(booleanStep3.y, and1.u[2]) annotation(
        Line(points = {{277, -173}, {285, -173}, {285, -155}, {288, -155}}, color = {255, 0, 255}));
      connect(booleanStep4.y, and1.u[1]) annotation(
        Line(points = {{279, -140}, {285, -140}, {285, -155}, {288, -155}}, color = {255, 0, 255}));
      connect(switch2.n, resistor3.p) annotation(
        Line(points = {{334, 21}, {334, 26}, {289, 26}, {289, 39}}, color = {0, 0, 255}));
      connect(resistor6.p, switch2.p) annotation(
        Line(points = {{342, 1}, {334, 1}}, color = {0, 0, 255}));
      connect(ground11.p, resistor6.n) annotation(
        Line(points = {{366, 2}, {362, 2}, {362, 1}}, color = {0, 0, 255}));
      connect(idealClosingSwitch.n, resistor4.p) annotation(
        Line(points = {{337, -37}, {288, -37}, {288, -21}}, color = {0, 0, 255}));
      connect(idealClosingSwitch.p, resistor7.p) annotation(
        Line(points = {{337, -57}, {344, -57}}, color = {0, 0, 255}));
      connect(resistor7.n, ground12.p) annotation(
        Line(points = {{364, -57}, {368, -57}, {368, -56}}, color = {0, 0, 255}));
      connect(idealClosingSwitch1.n, resistor5.p) annotation(
        Line(points = {{337, -95}, {284, -95}, {284, -81}, {288, -81}}, color = {0, 0, 255}));
      connect(idealClosingSwitch1.p, resistor8.p) annotation(
        Line(points = {{337, -115}, {344, -115}, {344, -117}}, color = {0, 0, 255}));
      connect(resistor8.n, ground13.p) annotation(
        Line(points = {{364, -117}, {368, -117}, {368, -114}}, color = {0, 0, 255}));
      connect(and1.y, switch2.control) annotation(
        Line(points = {{301, -155}, {317, -155}, {317, 11}, {322, 11}}, color = {255, 0, 255}));
      connect(booleanStep.y, and2.u[1]) annotation(
        Line(points = {{205, -140}, {211, -140}, {211, -155}, {214, -155}}, color = {255, 0, 255}));
      connect(booleanStep1.y, and2.u[2]) annotation(
        Line(points = {{203, -173}, {211, -173}, {211, -155}, {214, -155}}, color = {255, 0, 255}));
      connect(and2.y, idealClosingSwitch.control) annotation(
        Line(points = {{227, -155}, {237, -155}, {237, -47}, {325, -47}}, color = {255, 0, 255}));
      connect(booleanStep2.y, and3.u[1]) annotation(
        Line(points = {{357, -157}, {363, -157}, {363, -172}, {366, -172}}, color = {255, 0, 255}));
      connect(booleanStep5.y, and3.u[2]) annotation(
        Line(points = {{355, -190}, {363, -190}, {363, -172}, {366, -172}}, color = {255, 0, 255}));
      connect(and3.y, idealClosingSwitch1.control) annotation(
        Line(points = {{379, -172}, {384, -172}, {384, -135}, {320, -135}, {320, -105}, {325, -105}}, color = {255, 0, 255}));
      connect(constant16.y, ib_pu.u2) annotation(
        Line(points = {{250, 241}, {258, 241}, {258, 282}, {268, 282}}, color = {0, 0, 127}));
      connect(constant16.y, ia_pu.u2) annotation(
        Line(points = {{250, 241}, {258, 241}, {258, 309}, {268, 309}}, color = {0, 0, 127}));
      connect(realExpression45.y, ic_pu.u1) annotation(
        Line(points = {{230, 272}, {261, 272}, {261, 264}, {268, 264}}, color = {0, 0, 127}));
      connect(realExpression46.y, ia_pu.u1) annotation(
        Line(points = {{230, 304}, {245, 304}, {245, 317}, {268, 317}}, color = {0, 0, 127}));
      connect(constant16.y, ic_pu.u2) annotation(
        Line(points = {{250, 241}, {258, 241}, {258, 256}, {268, 256}}, color = {0, 0, 127}));
      connect(realExpression44.y, ib_pu.u1) annotation(
        Line(points = {{230, 288}, {268, 288}, {268, 290}}, color = {0, 0, 127}));
  connect(inversePark_fz.alpha, clarkInv.A) annotation(
        Line(points = {{485, 145}, {517, 145}, {517, 146}}, color = {0, 0, 127}));
  connect(inversePark_fz.beta, clarkInv.B) annotation(
        Line(points = {{485, 137}, {517, 137}, {517, 140}}, color = {0, 0, 127}));
  connect(constant2.y, clarkInv.C) annotation(
        Line(points = {{500, 98}, {517, 98}, {517, 135}}, color = {0, 0, 127}));
      annotation(
        Diagram(coordinateSystem(extent = {{-150, 300}, {700, -130}}, grid = {1, 1})),
        Icon(coordinateSystem(extent = {{-120, -120}, {480, 160}}, grid = {1, 1})));
    end ART52;
    
    model ClarkInv
      Modelica.Blocks.Interfaces.RealInput A annotation(
        Placement(transformation(extent = {{-148, 30}, {-108, 70}})));
      Modelica.Blocks.Interfaces.RealInput B annotation(
        Placement(transformation(extent = {{-150, -26}, {-110, 14}})));
      Modelica.Blocks.Interfaces.RealInput C annotation(
        Placement(transformation(extent = {{-148, -78}, {-108, -38}})));
      Modelica.Blocks.Interfaces.RealOutput a annotation(
        Placement(transformation(extent = {{100, 38}, {120, 58}})));
      Modelica.Blocks.Interfaces.RealOutput b annotation(
        Placement(transformation(extent = {{100, -12}, {120, 8}})));
      Modelica.Blocks.Interfaces.RealOutput c annotation(
        Placement(transformation(extent = {{100, -58}, {120, -38}})));
    equation
      a = A + C;
      b = (1/2)*(-A + sqrt(3)*B + 2*C);
      c = (1/2)*(-A - sqrt(3)*B + C);
      annotation(
        Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}})));
    end ClarkInv;
  end Testes;
  annotation(
    uses(iPSL(version = "0.8.1"), Modelica(version = "3.2.2"), PVSystems(version = "0.6.3")),
    Documentation(info = "<html>
<table cellspacing=\"2\" cellpadding=\"0\" border=\"0\"<tr>
<td><p>Reference</p></td>
<td><p>EMTP-RV</p></td>
</tr>
<tr>
<td><p>Last update</p></td>
<td><p>30/04/2017</p></td>
</tr>
<tr>
<td><p>Author</p></td>
<td><p>MAA Murad and Luigi Vanfretti, SmarTS Lab, KTH Royal Institute of Technology</p></td>
</tr>
<tr>
<td><p>Contact</p></td>
<td><p><a href=\"mailto:luigiv@kth.se\">luigiv@kth.se</a></p></td>
</tr>
</table>
</html>"));
end VSC_FZ;
