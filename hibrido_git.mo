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
        Placement(visible = true, transformation(origin = {0, -52}, extent = {{100, -10}, {120, 10}}, rotation = 0), iconTransformation(origin = {0, -30}, extent = {{100, -10}, {120, 10}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealOutput q annotation(
        Placement(visible = true, transformation(origin = {0, -82}, extent = {{100, -10}, {120, 10}}, rotation = 0), iconTransformation(origin = {0, -90}, extent = {{100, -10}, {120, 10}}, rotation = 0)));
      VSC_FZ.Testes.PI1 pi1(KI = 87.73, KP = 0.9872, init_value = 0, max = 1000, min = -1000) annotation(
        Placement(visible = true, transformation(origin = {6, 24}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealOutput freqout annotation(
        Placement(visible = true, transformation(origin = {0, 50}, extent = {{100, -10}, {120, 10}}, rotation = 0), iconTransformation(origin = {0, 22}, extent = {{100, -10}, {120, 10}}, rotation = 0)));
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
      connect(add.y, freqout) annotation(
        Line(points = {{52, 58}, {94, 58}, {94, 50}, {110, 50}}, color = {0, 0, 127}));
      annotation(
        Diagram(coordinateSystem(extent = {{-140, 80}, {120, -100}})),
        Icon(graphics = {Line(points = {{-70, 0}, {-50, 60}, {-30, 0}, {-10, -60}, {10, 0}, {30, 60}, {50, 0}, {70, -60}, {90, 0}}, color = {0, 0, 255}, smooth = Smooth.Bezier), Line(points = {{-90, 0}, {-64, 60}, {-44, 0}, {-18, -60}, {2, 0}, {22, 60}, {44, 0}, {64, -60}, {88, 0}}, color = {255, 0, 0}, smooth = Smooth.Bezier)}, coordinateSystem(extent = {{-100, -100}, {100, 100}})),
        Documentation(info = "<html>
                <p>
                  Phase-locked loop. Given a sinusoidal input, extract the phase.
                </p>
              </html>"));
    end PLL_fz2;

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

    model D_D_Transformer
      parameter Real R_trans_prim;
      parameter Real L_trans_prim;
      parameter Real R_trans_secon;
      parameter Real L_trans_secon;
      parameter Real Vac_primary;
      parameter Real Vac_secondary;
      parameter Real trans_ratio;
      parameter Real Lm1;
      parameter Boolean considerMagnetization;
      Modelica.Electrical.Analog.Ideal.IdealTransformer transformer(Lm1 = Lm1, considerMagnetization = considerMagnetization, n = trans_ratio) annotation(
        Placement(visible = true, transformation(origin = {-2, 54}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Resistor resistor(R = R_trans_prim) annotation(
        Placement(visible = true, transformation(origin = {-62, 64}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Inductor inductor(L = L_trans_prim) annotation(
        Placement(visible = true, transformation(origin = {-32, 64}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Resistor resistor1(R = R_trans_secon) annotation(
        Placement(transformation(origin = {32, 64}, extent = {{-10, -10}, {10, 10}})));
      Modelica.Electrical.Analog.Basic.Inductor inductor1(L = L_trans_secon) annotation(
        Placement(visible = true, transformation(origin = {62, 64}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Inductor inductor2(L = L_trans_secon) annotation(
        Placement(visible = true, transformation(origin = {64, 12}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Resistor resistor2(R = R_trans_secon) annotation(
        Placement(visible = true, transformation(origin = {34, 12}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Resistor resistor3(R = R_trans_prim) annotation(
        Placement(visible = true, transformation(origin = {-60, 12}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Ideal.IdealTransformer idealTransformer(Lm1 = Lm1, considerMagnetization = considerMagnetization, n = trans_ratio) annotation(
        Placement(visible = true, transformation(origin = {0, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Inductor inductor3(L = L_trans_prim) annotation(
        Placement(visible = true, transformation(origin = {-30, 12}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Inductor inductor4(L = L_trans_secon) annotation(
        Placement(visible = true, transformation(origin = {64, -38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Resistor resistor4(R = R_trans_secon) annotation(
        Placement(visible = true, transformation(origin = {34, -38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Resistor resistor5(R = R_trans_prim) annotation(
        Placement(visible = true, transformation(origin = {-60, -38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Ideal.IdealTransformer idealTransformer1(Lm1 = Lm1, considerMagnetization = considerMagnetization, n = trans_ratio) annotation(
        Placement(visible = true, transformation(origin = {0, -48}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Inductor inductor5(L = L_trans_prim) annotation(
        Placement(visible = true, transformation(origin = {-30, -38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Interfaces.Pin pin annotation(
        Placement(visible = true, transformation(origin = {-106, 58}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-106, 58}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Interfaces.Pin pin1 annotation(
        Placement(visible = true, transformation(origin = {-106, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-106, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Interfaces.Pin pin2 annotation(
        Placement(visible = true, transformation(origin = {-104, -38}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-104, -38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Interfaces.Pin pin3 annotation(
        Placement(visible = true, transformation(origin = {106, 62}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {106, 62}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Interfaces.Pin pin4 annotation(
        Placement(visible = true, transformation(origin = {106, 12}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {106, 12}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Interfaces.Pin pin5 annotation(
        Placement(visible = true, transformation(origin = {104, -38}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {104, -38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.VoltageSensor v1d annotation(
        Placement(visible = true, transformation(origin = {-28, 35}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.VoltageSensor vin annotation(
        Placement(visible = true, transformation(origin = {-96, 37}, extent = {{-7, -7}, {7, 7}}, rotation = -90)));
      Modelica.Electrical.Analog.Sensors.VoltageSensor vout annotation(
        Placement(visible = true, transformation(origin = {96, 39}, extent = {{7, -7}, {-7, 7}}, rotation = 90)));
      Modelica.Electrical.Analog.Sensors.VoltageSensor voltageSensor annotation(
        Placement(visible = true, transformation(origin = {-8, 93}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Ground ground annotation(
        Placement(visible = true, transformation(origin = {8, 80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.VoltageSensor voltageSensor1 annotation(
        Placement(visible = true, transformation(origin = {40, 93}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Ground ground1 annotation(
        Placement(visible = true, transformation(origin = {56, 80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(resistor.n, inductor.p) annotation(
        Line(points = {{-52, 64}, {-42, 64}}, color = {0, 0, 255}));
      connect(inductor.n, transformer.p1) annotation(
        Line(points = {{-22, 64}, {-12, 64}}, color = {0, 0, 255}));
      connect(transformer.p2, resistor1.p) annotation(
        Line(points = {{8, 64}, {22, 64}}, color = {0, 0, 255}));
      connect(resistor1.n, inductor1.p) annotation(
        Line(points = {{42, 64}, {52, 64}}, color = {0, 0, 255}));
      connect(resistor2.n, inductor2.p) annotation(
        Line(points = {{44, 12}, {54, 12}}, color = {0, 0, 255}));
      connect(inductor3.n, idealTransformer.p1) annotation(
        Line(points = {{-20, 12}, {-10, 12}}, color = {0, 0, 255}));
      connect(idealTransformer.p2, resistor2.p) annotation(
        Line(points = {{10, 12}, {24, 12}}, color = {0, 0, 255}));
      connect(resistor3.n, inductor3.p) annotation(
        Line(points = {{-50, 12}, {-40, 12}}, color = {0, 0, 255}));
      connect(resistor4.n, inductor4.p) annotation(
        Line(points = {{44, -38}, {54, -38}}, color = {0, 0, 255}));
      connect(inductor5.n, idealTransformer1.p1) annotation(
        Line(points = {{-20, -38}, {-10, -38}}, color = {0, 0, 255}));
      connect(idealTransformer1.p2, resistor4.p) annotation(
        Line(points = {{10, -38}, {24, -38}}, color = {0, 0, 255}));
      connect(resistor5.n, inductor5.p) annotation(
        Line(points = {{-50, -38}, {-40, -38}}, color = {0, 0, 255}));
      connect(transformer.n1, resistor3.p) annotation(
        Line(points = {{-12, 44}, {-70, 44}, {-70, 12}}, color = {0, 0, 255}));
      connect(idealTransformer.n1, resistor5.p) annotation(
        Line(points = {{-10, -8}, {-70, -8}, {-70, -38}}, color = {0, 0, 255}));
      connect(idealTransformer1.n1, resistor.p) annotation(
        Line(points = {{-10, -58}, {-80, -58}, {-80, 64}, {-72, 64}}, color = {0, 0, 255}));
      connect(transformer.n2, inductor2.n) annotation(
        Line(points = {{8, 44}, {74, 44}, {74, 12}}, color = {0, 0, 255}));
      connect(idealTransformer.n2, inductor4.n) annotation(
        Line(points = {{10, -8}, {74, -8}, {74, -38}}, color = {0, 0, 255}));
      connect(idealTransformer1.n2, inductor1.n) annotation(
        Line(points = {{10, -58}, {84, -58}, {84, 64}, {72, 64}}, color = {0, 0, 255}));
      connect(inductor1.n, pin3) annotation(
        Line(points = {{72, 64}, {106, 64}, {106, 62}}, color = {0, 0, 255}));
      connect(inductor2.n, pin4) annotation(
        Line(points = {{74, 12}, {106, 12}}, color = {0, 0, 255}));
      connect(inductor4.n, pin5) annotation(
        Line(points = {{74, -38}, {104, -38}}, color = {0, 0, 255}));
      connect(resistor.p, pin) annotation(
        Line(points = {{-72, 64}, {-94, 64}, {-94, 58}, {-106, 58}}, color = {0, 0, 255}));
      connect(resistor3.p, pin1) annotation(
        Line(points = {{-70, 12}, {-106, 12}, {-106, 10}}, color = {0, 0, 255}));
      connect(resistor5.p, pin2) annotation(
        Line(points = {{-70, -38}, {-104, -38}}, color = {0, 0, 255}));
      connect(v1d.p, transformer.p1) annotation(
        Line(points = {{-34, 36}, {-42, 36}, {-42, 50}, {-16, 50}, {-16, 57}, {-12, 57}, {-12, 64}}, color = {0, 0, 255}));
      connect(v1d.n, transformer.n1) annotation(
        Line(points = {{-20, 36}, {-12, 36}, {-12, 44}}, color = {0, 0, 255}));
      connect(vin.p, resistor.p) annotation(
        Line(points = {{-96, 44}, {-94, 44}, {-94, 64}, {-72, 64}}, color = {0, 0, 255}));
      connect(vin.n, pin1) annotation(
        Line(points = {{-96, 30}, {-96, 10}, {-106, 10}}, color = {0, 0, 255}));
      connect(vout.p, pin3) annotation(
        Line(points = {{96, 46}, {96, 62}, {106, 62}}, color = {0, 0, 255}));
      connect(vout.n, pin4) annotation(
        Line(points = {{96, 32}, {96, 12}, {106, 12}}, color = {0, 0, 255}));
      connect(voltageSensor.p, transformer.p1) annotation(
        Line(points = {{-14, 94}, {-20, 94}, {-20, 70}, {-12, 70}, {-12, 64}}, color = {0, 0, 255}));
      connect(voltageSensor.n, ground.p) annotation(
        Line(points = {{0, 94}, {8, 94}, {8, 90}}, color = {0, 0, 255}));
      connect(voltageSensor1.n, ground1.p) annotation(
        Line(points = {{47, 93}, {55, 93}, {55, 89}}, color = {0, 0, 255}));
      connect(voltageSensor1.p, resistor1.p) annotation(
        Line(points = {{34, 94}, {34, 96}, {22, 96}, {22, 64}}, color = {0, 0, 255}));
    end D_D_Transformer;

    model transformer_error
      Modelica.Blocks.Interfaces.RealInput p1 annotation(
        Placement(visible = true, transformation(origin = {-120, 32}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-106, 36}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealInput p2 annotation(
        Placement(visible = true, transformation(origin = {-120, -30}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-110, -32}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealOutput error annotation(
        Placement(visible = true, transformation(origin = {108, 4}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {106, 4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      error = p1 - p2;
      annotation(
        uses(Modelica(version = "3.2.3")));
    end transformer_error;

    model mod
      Modelica.Blocks.Interfaces.RealInput I annotation(
        Placement(visible = true, transformation(origin = {-117, 59}, extent = {{-17, -17}, {17, 17}}, rotation = 0), iconTransformation(origin = {-112, 46}, extent = {{-12, -12}, {12, 12}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealOutput O annotation(
        Placement(visible = true, transformation(origin = {0, 2}, extent = {{100, -10}, {120, 10}}, rotation = 0), iconTransformation(origin = {0, 2}, extent = {{100, -10}, {120, 10}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealInput Ref annotation(
        Placement(visible = true, transformation(origin = {-117, -41}, extent = {{-17, -17}, {17, 17}}, rotation = 0), iconTransformation(origin = {-112, -46}, extent = {{-12, -12}, {12, 12}}, rotation = 0)));
    equation
      O = I - Ref*floor(I/Ref);
    end mod;

    model SOGI
      parameter Real frequency;
      parameter Real k;
      Modelica.Blocks.Interfaces.RealOutput sigout annotation(
        Placement(visible = true, transformation(origin = {30, 34}, extent = {{100, -10}, {120, 10}}, rotation = 0), iconTransformation(origin = {0, 40}, extent = {{100, -10}, {120, 10}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealInput sigin annotation(
        Placement(visible = true, transformation(origin = {-139, 39}, extent = {{-17, -17}, {17, 17}}, rotation = 0), iconTransformation(origin = {-112, 2.22045e-16}, extent = {{-12, -12}, {12, 12}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealOutput qsigout annotation(
        Placement(visible = true, transformation(origin = {0, -76}, extent = {{100, -10}, {120, 10}}, rotation = 0), iconTransformation(origin = {0, -40}, extent = {{100, -10}, {120, 10}}, rotation = 0)));
      Modelica.Blocks.Continuous.Integrator integrator annotation(
        Placement(visible = true, transformation(origin = {28, 34}, extent = {{60, -10}, {80, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Add add1(k1 = -1) annotation(
        Placement(visible = true, transformation(origin = {-116, -12}, extent = {{30, 48}, {50, 68}}, rotation = 0)));
      Modelica.Blocks.Math.Gain gain annotation(
        Placement(visible = true, transformation(origin = {-32, 46}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Add add2(k1 = +1, k2 = -1) annotation(
        Placement(visible = true, transformation(origin = {-34, -18}, extent = {{30, 48}, {50, 68}}, rotation = 0)));
      Modelica.Blocks.Math.Product product annotation(
        Placement(visible = true, transformation(origin = {50, 34}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Product product1 annotation(
        Placement(visible = true, transformation(origin = {48, -40}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      Modelica.Blocks.Continuous.Integrator integrator1 annotation(
        Placement(visible = true, transformation(origin = {168, -46}, extent = {{-60, -10}, {-80, 10}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealInput freqin annotation(
        Placement(visible = true, transformation(origin = {-139, 1}, extent = {{-17, -17}, {17, 17}}, rotation = 0), iconTransformation(origin = {-112, 60}, extent = {{-12, -12}, {12, 12}}, rotation = 0)));
    equation
      connect(integrator.y, sigout) annotation(
        Line(points = {{109, 34}, {140, 34}}, color = {0, 0, 127}));
      connect(sigin, add1.u2) annotation(
        Line(points = {{-139, 39}, {-89, 39}}, color = {0, 0, 127}));
      connect(add1.y, gain.u) annotation(
        Line(points = {{-65, 46}, {-45, 46}}, color = {0, 0, 127}));
      connect(gain.y, add2.u1) annotation(
        Line(points = {{-21, 46}, {-7, 46}}, color = {0, 0, 127}));
      connect(add2.y, product.u1) annotation(
        Line(points = {{17, 40}, {37, 40}}, color = {0, 0, 127}));
      connect(product.y, integrator.u) annotation(
        Line(points = {{62, 34}, {86, 34}}, color = {0, 0, 127}));
      connect(integrator1.y, product1.u2) annotation(
        Line(points = {{87, -46}, {60, -46}}, color = {0, 0, 127}));
      connect(product1.y, qsigout) annotation(
        Line(points = {{37, -40}, {23, -40}, {23, -76}, {109, -76}}, color = {0, 0, 127}));
      connect(product1.y, add2.u2) annotation(
        Line(points = {{37, -40}, {-34, -40}, {-34, 26}, {-16, 26}, {-16, 34}, {-6, 34}}, color = {0, 0, 127}));
      connect(integrator.y, integrator1.u) annotation(
        Line(points = {{110, 34}, {120, 34}, {120, -46}, {110, -46}}, color = {0, 0, 127}));
      connect(integrator.y, add1.u1) annotation(
        Line(points = {{110, 34}, {120, 34}, {120, 72}, {-100, 72}, {-100, 52}, {-88, 52}}, color = {0, 0, 127}));
      connect(freqin, product.u2) annotation(
        Line(points = {{-139, 1}, {28, 1}, {28, 28}, {38, 28}}, color = {0, 0, 127}));
      connect(freqin, product1.u1) annotation(
        Line(points = {{-139, 1}, {66, 1}, {66, -34}, {60, -34}}, color = {0, 0, 127}));
    end SOGI;

    model DSOGI_PLL
      VSC_FZ.Testes.Clark3p clark3p annotation(
        Placement(visible = true, transformation(origin = {-98, 32}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      VSC_FZ.Testes.SOGI sogi(frequency = 60, k = 0.7) annotation(
        Placement(visible = true, transformation(origin = {-54, 44}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      VSC_FZ.Testes.SOGI sogi1(frequency = 60, k = 0.7) annotation(
        Placement(visible = true, transformation(origin = {-54, 14}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Add add1(k1 = +1, k2 = -1) annotation(
        Placement(visible = true, transformation(origin = {-48, -6}, extent = {{30, 48}, {50, 68}}, rotation = 0)));
      Modelica.Blocks.Math.Gain gain(k = 0.5) annotation(
        Placement(visible = true, transformation(origin = {38, 52}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Add add(k1 = +1, k2 = +1) annotation(
        Placement(visible = true, transformation(origin = {-48, -34}, extent = {{30, 48}, {50, 68}}, rotation = 0)));
      Modelica.Blocks.Math.Gain gain1(k = 0.5) annotation(
        Placement(visible = true, transformation(origin = {38, 24}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      VSC_FZ.Testes.PLL_fz2 pLL_fz2(frequency = 60.08)  annotation(
        Placement(visible = true, transformation(origin = {84, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      VSC_FZ.Testes.mod mod annotation(
        Placement(visible = true, transformation(origin = {128, 4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant lineFreq(k = 2*Modelica.Constants.pi) annotation(
        Placement(visible = true, transformation(origin = {74, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealInput V_abc[3] annotation(
        Placement(visible = true, transformation(origin = {-28, 62}, extent = {{-146, -50}, {-106, -10}}, rotation = 0), iconTransformation(origin = {-111, 1}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealOutput d annotation(
        Placement(visible = true, transformation(origin = {62, 50}, extent = {{100, -10}, {120, 10}}, rotation = 0), iconTransformation(origin = {0, 90}, extent = {{100, -10}, {120, 10}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealOutput q annotation(
        Placement(visible = true, transformation(origin = {62, 30}, extent = {{100, -10}, {120, 10}}, rotation = 0), iconTransformation(origin = {0, 30}, extent = {{100, -10}, {120, 10}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealOutput ang annotation(
        Placement(visible = true, transformation(origin = {60, -12}, extent = {{100, -10}, {120, 10}}, rotation = 0), iconTransformation(origin = {0, -90}, extent = {{100, -10}, {120, 10}}, rotation = 0)));
    equation
      connect(clark3p.V_ab_D[1], sogi.sigin) annotation(
        Line(points = {{-87, 32}, {-77, 32}, {-77, 44}, {-67, 44}}, color = {0, 0, 127}));
      connect(sogi1.sigin, clark3p.V_ab_D[2]) annotation(
        Line(points = {{-65.2, 14}, {-75.2, 14}, {-75.2, 32}, {-85.2, 32}}, color = {0, 0, 127}));
      connect(sogi.sigout, add1.u1) annotation(
        Line(points = {{-43, 48}, {-33, 48}, {-33, 58}, {-21, 58}}, color = {0, 0, 127}));
      connect(sogi.qsigout, add.u1) annotation(
        Line(points = {{-43, 40}, {-31, 40}, {-31, 30}, {-21, 30}}, color = {0, 0, 127}));
      connect(sogi1.sigout, add.u2) annotation(
        Line(points = {{-43, 18}, {-21, 18}}, color = {0, 0, 127}));
      connect(add1.u2, sogi1.qsigout) annotation(
        Line(points = {{-20, 46}, {-26, 46}, {-26, 10}, {-42, 10}}, color = {0, 0, 127}));
      connect(add1.y, gain.u) annotation(
        Line(points = {{3, 52}, {25, 52}}, color = {0, 0, 127}));
      connect(add.y, gain1.u) annotation(
        Line(points = {{3, 24}, {25, 24}}, color = {0, 0, 127}));
      connect(gain.y, pLL_fz2.alpha) annotation(
        Line(points = {{49, 52}, {61, 52}, {61, 42}, {71, 42}}, color = {0, 0, 127}));
      connect(gain1.y, pLL_fz2.beta) annotation(
        Line(points = {{49, 24}, {63, 24}, {63, 34}, {71, 34}}, color = {0, 0, 127}));
      connect(pLL_fz2.freqout, sogi.freqin) annotation(
        Line(points = {{95, 45}, {103, 45}, {103, 71}, {-75, 71}, {-75, 49}, {-67, 49}}, color = {0, 0, 127}));
      connect(sogi1.freqin, sogi.freqin) annotation(
        Line(points = {{-65.2, 20}, {-73.2, 20}, {-73.2, 50}, {-65.2, 50}}, color = {0, 0, 127}));
      connect(lineFreq.y, mod.Ref) annotation(
        Line(points = {{85, 0}, {116, 0}, {116, -1}, {117, -1}}, color = {0, 0, 127}));
      connect(V_abc[1], clark3p.V_abc_D[1]) annotation(
        Line(points = {{-154, 32}, {-110, 32}}, color = {0, 0, 127}, thickness = 0.5));
      connect(V_abc[2], clark3p.V_abc_D[2]) annotation(
        Line(points = {{-154, 32}, {-110, 32}}, color = {0, 0, 127}, thickness = 0.5));
      connect(V_abc[3], clark3p.V_abc_D[3]) annotation(
        Line(points = {{-154, 32}, {-110, 32}}, color = {0, 0, 127}, thickness = 0.5));
      connect(pLL_fz2.d, d) annotation(
        Line(points = {{96, 36}, {144, 36}, {144, 50}, {172, 50}}, color = {0, 0, 127}));
      connect(pLL_fz2.q, q) annotation(
        Line(points = {{96, 30}, {172, 30}}, color = {0, 0, 127}));
      connect(pLL_fz2.theta, mod.I) annotation(
        Line(points = {{96, 46}, {108, 46}, {108, 8}, {116, 8}}, color = {0, 0, 127}));
      connect(pLL_fz2.theta, ang) annotation(
        Line(points = {{96, 46}, {148, 46}, {148, -12}, {170, -12}}, color = {0, 0, 127}));
      annotation(
        Diagram(coordinateSystem(extent = {{-180, 80}, {180, -20}})));
    end DSOGI_PLL;

    model Res
      parameter Real freq;
      parameter Real h;
      parameter Real KIh;
      Modelica.Blocks.Interfaces.RealOutput sigout annotation(
        Placement(visible = true, transformation(origin = {0, 40}, extent = {{100, -10}, {120, 10}}, rotation = 0), iconTransformation(origin = {0, 40}, extent = {{100, -10}, {120, 10}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealInput sigin annotation(
        Placement(visible = true, transformation(origin = {-100, 46}, extent = {{-12, -12}, {12, 12}}, rotation = 0), iconTransformation(origin = {-112, 0}, extent = {{-12, -12}, {12, 12}}, rotation = 0)));
      Modelica.Blocks.Continuous.Integrator integrator annotation(
        Placement(visible = true, transformation(origin = {-22, 40}, extent = {{60, -10}, {80, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Gain gain(k = KIh) annotation(
        Placement(visible = true, transformation(origin = {-32, 46}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Add add2(k1 = +1, k2 = -1) annotation(
        Placement(visible = true, transformation(origin = {-34, -18}, extent = {{30, 48}, {50, 68}}, rotation = 0)));
      Modelica.Blocks.Math.Product product1 annotation(
        Placement(visible = true, transformation(origin = {10, -4}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      Modelica.Blocks.Continuous.Integrator integrator1 annotation(
        Placement(visible = true, transformation(origin = {144, 2}, extent = {{-60, -10}, {-80, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant const(k = h*h*freq*freq) annotation(
        Placement(visible = true, transformation(origin = {10, -38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(gain.y, add2.u1) annotation(
        Line(points = {{-21, 46}, {-7, 46}}, color = {0, 0, 127}));
      connect(product1.y, add2.u2) annotation(
        Line(points = {{-1, -4}, {-16, -4}, {-16, 34}, {-6, 34}}, color = {0, 0, 127}));
      connect(integrator.y, integrator1.u) annotation(
        Line(points = {{59, 40}, {100, 40}, {100, 2}, {86, 2}}, color = {0, 0, 127}));
      connect(add2.y, integrator.u) annotation(
        Line(points = {{18, 40}, {36, 40}}, color = {0, 0, 127}));
      connect(product1.u1, integrator1.y) annotation(
        Line(points = {{22, 2}, {64, 2}}, color = {0, 0, 127}));
      connect(product1.u2, const.y) annotation(
        Line(points = {{22, -10}, {32, -10}, {32, -38}, {21, -38}}, color = {0, 0, 127}));
      connect(sigout, integrator.y) annotation(
        Line(points = {{110, 40}, {60, 40}}, color = {0, 0, 127}));
      connect(sigin, gain.u) annotation(
        Line(points = {{-100, 46}, {-44, 46}}, color = {0, 0, 127}));
      annotation(
        Diagram(coordinateSystem(extent = {{-120, 60}, {120, -60}})));
    end Res;

    model InnerControl_PIR
      parameter Real Ictrl_KP;
      parameter Real Ictrl_KI;
      parameter Real Lac_eq_pu;
      parameter Real initvd;
      parameter Real initvq;
      parameter Real r;
      parameter Real lw;
      parameter Real freq;
      parameter Real h;
      parameter Real KIh;
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
        Placement(visible = true, transformation(origin = {14, 8}, extent = {{18, 46}, {38, 66}}, rotation = 0)));
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
      VSC_FZ.Testes.PI1 pi1(KI = Ictrl_KI, KP = Ictrl_KP, init_value = 0, max = 10000, min = -10000, time_step = 0.00001) annotation(
        Placement(visible = true, transformation(origin = {-22, 74}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      VSC_FZ.Testes.PI1 pi11(KI = Ictrl_KI, KP = Ictrl_KP, init_value = 0, max = 10000, min = -10000, time_step = 0.00001) annotation(
        Placement(visible = true, transformation(origin = {-26, -46}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Add add2 annotation(
        Placement(visible = true, transformation(origin = {80, 58}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
      Modelica.Blocks.Math.Product product2 annotation(
        Placement(visible = true, transformation(origin = {84, 8}, extent = {{-46, 14}, {-34, 26}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant1(k = r) annotation(
        Placement(visible = true, transformation(origin = {100, 10}, extent = {{-88, -8}, {-74, 6}}, rotation = 0)));
      Modelica.Blocks.Math.Product product3 annotation(
        Placement(visible = true, transformation(origin = {82, -36}, extent = {{-46, 14}, {-34, 26}}, rotation = 0)));
      Modelica.Blocks.Math.Add add3 annotation(
        Placement(visible = true, transformation(origin = {74, -48}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
      Res res(KIh = KIh, freq = freq, h = h) annotation(
        Placement(visible = true, transformation(origin = {-20, 48}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Add add4 annotation(
        Placement(visible = true, transformation(origin = {8, 64}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
      Res res1(KIh = KIh, freq = freq, h = h) annotation(
        Placement(visible = true, transformation(origin = {-24, -70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Add add5 annotation(
        Placement(visible = true, transformation(origin = {0, -52}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
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
      connect(product.y, add3_2.u1) annotation(
        Line(points = {{-33.4, 20}, {-6, 20}, {-6, -38}, {18, -38}}, color = {0, 0, 127}));
      connect(add.y, pi1.u) annotation(
        Line(points = {{-52, 60}, {-43, 60}, {-43, 74}, {-34, 74}}, color = {0, 0, 127}));
      connect(add1.y, pi11.u) annotation(
        Line(points = {{-56, -46}, {-38, -46}}, color = {0, 0, 127}));
      connect(product2.u1, Id) annotation(
        Line(points = {{37, 32}, {-90, 32}, {-90, 18}, {-124, 18}}, color = {0, 0, 127}));
      connect(add3_1.y, add2.u1) annotation(
        Line(points = {{53, 64}, {54, 64}, {54, 63}, {70, 63}}, color = {0, 0, 127}));
      connect(product2.y, add2.u2) annotation(
        Line(points = {{50, 28}, {54, 28}, {54, 53}, {70, 53}}, color = {0, 0, 127}));
      connect(product2.u2, constant1.y) annotation(
        Line(points = {{36, 24}, {32, 24}, {32, 9}, {27, 9}}, color = {0, 0, 127}));
      connect(add2.y, vdref) annotation(
        Line(points = {{89, 58}, {97.5, 58}, {97.5, 56}, {110, 56}}, color = {0, 0, 127}));
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
      connect(add3_2.u3, vq) annotation(
        Line(points = {{18, -54}, {14, -54}, {14, -94}, {-126, -94}}, color = {0, 0, 127}));
      connect(product1.y, add3_1.u3) annotation(
        Line(points = {{-28, -18}, {0, -18}, {0, 44}, {24, 44}, {24, 56}, {30, 56}}, color = {0, 0, 127}));
      connect(add3_1.u1, vd) annotation(
        Line(points = {{30, 72}, {22, 72}, {22, 90}, {-122, 90}}, color = {0, 0, 127}));
      connect(add.y, res.sigin) annotation(
        Line(points = {{-52, 60}, {-44, 60}, {-44, 48}, {-32, 48}}, color = {0, 0, 127}));
      connect(pi1.y, add4.u1) annotation(
        Line(points = {{-10, 74}, {-6, 74}, {-6, 68}, {-2, 68}}, color = {0, 0, 127}));
      connect(res.sigout, add4.u2) annotation(
        Line(points = {{-8, 52}, {-8, 60}, {-2, 60}}, color = {0, 0, 127}));
      connect(add4.y, add3_1.u2) annotation(
        Line(points = {{16, 64}, {30, 64}}, color = {0, 0, 127}));
      connect(pi11.y, add5.u1) annotation(
        Line(points = {{-14, -46}, {-10, -46}, {-10, -48}}, color = {0, 0, 127}));
      connect(res1.sigout, add5.u2) annotation(
        Line(points = {{-12, -66}, {-12, -56}, {-10, -56}}, color = {0, 0, 127}));
      connect(add5.y, add3_2.u2) annotation(
        Line(points = {{8, -52}, {10, -52}, {10, -46}, {18, -46}}, color = {0, 0, 127}));
      connect(res1.sigin, add1.y) annotation(
        Line(points = {{-36, -70}, {-46, -70}, {-46, -46}, {-56, -46}}, color = {0, 0, 127}));
      annotation(
        Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}})));
    end InnerControl_PIR;

    model Yn_d_Transformer
      parameter Real R_trans_prim;
      parameter Real L_trans_prim;
      parameter Real R_trans_secon;
      parameter Real L_trans_secon;
      parameter Real Vac_primary;
      parameter Real Vac_secondary;
      parameter Real trans_ratio = (Vac_primary/Vac_secondary);
      parameter Real Lmag;
      parameter Real Rcor;
      Modelica.Electrical.Analog.Ideal.IdealTransformer transformer(Lm1 = 0, considerMagnetization = false, n = trans_ratio) annotation(
        Placement(visible = true, transformation(origin = {-2, 54}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Resistor resistor(R = R_trans_prim) annotation(
        Placement(visible = true, transformation(origin = {-62, 64}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Inductor inductor(L = L_trans_prim) annotation(
        Placement(visible = true, transformation(origin = {-32, 64}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Resistor resistor1(R = R_trans_secon) annotation(
        Placement(transformation(origin = {32, 64}, extent = {{-10, -10}, {10, 10}})));
      Modelica.Electrical.Analog.Basic.Inductor inductor1(L = L_trans_secon) annotation(
        Placement(visible = true, transformation(origin = {62, 64}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Inductor inductor2(L = L_trans_secon) annotation(
        Placement(visible = true, transformation(origin = {64, 12}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Resistor resistor2(R = R_trans_secon) annotation(
        Placement(visible = true, transformation(origin = {34, 12}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Resistor resistor3(R = R_trans_prim) annotation(
        Placement(visible = true, transformation(origin = {-60, 12}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Ideal.IdealTransformer idealTransformer(Lm1 = 0, considerMagnetization = false, n = trans_ratio) annotation(
        Placement(visible = true, transformation(origin = {0, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Inductor inductor3(L = L_trans_prim) annotation(
        Placement(visible = true, transformation(origin = {-30, 12}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Inductor inductor4(L = L_trans_secon) annotation(
        Placement(visible = true, transformation(origin = {64, -38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Resistor resistor4(R = R_trans_secon) annotation(
        Placement(visible = true, transformation(origin = {34, -38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Resistor resistor5(R = R_trans_prim) annotation(
        Placement(visible = true, transformation(origin = {-60, -38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Ideal.IdealTransformer idealTransformer1(Lm1 = 0, considerMagnetization = false, n = trans_ratio) annotation(
        Placement(visible = true, transformation(origin = {0, -48}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Inductor inductor5(L = L_trans_prim) annotation(
        Placement(visible = true, transformation(origin = {-30, -38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Interfaces.Pin pin annotation(
        Placement(visible = true, transformation(origin = {-106, 58}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-106, 58}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Interfaces.Pin pin1 annotation(
        Placement(visible = true, transformation(origin = {-106, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-106, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Interfaces.Pin pin2 annotation(
        Placement(visible = true, transformation(origin = {-104, -38}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-104, -38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Interfaces.Pin pin3 annotation(
        Placement(visible = true, transformation(origin = {106, 62}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {106, 62}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Interfaces.Pin pin4 annotation(
        Placement(visible = true, transformation(origin = {106, 12}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {106, 12}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Interfaces.Pin pin5 annotation(
        Placement(visible = true, transformation(origin = {104, -38}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {104, -38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Interfaces.Pin pin6 annotation(
        Placement(visible = true, transformation(origin = {-104, -70}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-106, -76}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(resistor.n, inductor.p) annotation(
        Line(points = {{-52, 64}, {-42, 64}}, color = {0, 0, 255}));
      connect(inductor.n, transformer.p1) annotation(
        Line(points = {{-22, 64}, {-12, 64}}, color = {0, 0, 255}));
      connect(transformer.p2, resistor1.p) annotation(
        Line(points = {{8, 64}, {22, 64}}, color = {0, 0, 255}));
      connect(resistor1.n, inductor1.p) annotation(
        Line(points = {{42, 64}, {52, 64}}, color = {0, 0, 255}));
      connect(resistor2.n, inductor2.p) annotation(
        Line(points = {{44, 12}, {54, 12}}, color = {0, 0, 255}));
      connect(inductor3.n, idealTransformer.p1) annotation(
        Line(points = {{-20, 12}, {-10, 12}}, color = {0, 0, 255}));
      connect(idealTransformer.p2, resistor2.p) annotation(
        Line(points = {{10, 12}, {24, 12}}, color = {0, 0, 255}));
      connect(resistor3.n, inductor3.p) annotation(
        Line(points = {{-50, 12}, {-40, 12}}, color = {0, 0, 255}));
      connect(resistor4.n, inductor4.p) annotation(
        Line(points = {{44, -38}, {54, -38}}, color = {0, 0, 255}));
      connect(inductor5.n, idealTransformer1.p1) annotation(
        Line(points = {{-20, -38}, {-10, -38}}, color = {0, 0, 255}));
      connect(idealTransformer1.p2, resistor4.p) annotation(
        Line(points = {{10, -38}, {24, -38}}, color = {0, 0, 255}));
      connect(resistor5.n, inductor5.p) annotation(
        Line(points = {{-50, -38}, {-40, -38}}, color = {0, 0, 255}));
      connect(transformer.n2, inductor2.n) annotation(
        Line(points = {{8, 44}, {74, 44}, {74, 12}}, color = {0, 0, 255}));
      connect(idealTransformer.n2, inductor4.n) annotation(
        Line(points = {{10, -8}, {74, -8}, {74, -38}}, color = {0, 0, 255}));
      connect(idealTransformer1.n2, inductor1.n) annotation(
        Line(points = {{10, -58}, {84, -58}, {84, 64}, {72, 64}}, color = {0, 0, 255}));
      connect(inductor1.n, pin3) annotation(
        Line(points = {{72, 64}, {106, 64}, {106, 62}}, color = {0, 0, 255}));
      connect(inductor2.n, pin4) annotation(
        Line(points = {{74, 12}, {106, 12}}, color = {0, 0, 255}));
      connect(inductor4.n, pin5) annotation(
        Line(points = {{74, -38}, {104, -38}}, color = {0, 0, 255}));
      connect(resistor.p, pin) annotation(
        Line(points = {{-72, 64}, {-94, 64}, {-94, 58}, {-106, 58}}, color = {0, 0, 255}));
      connect(resistor5.p, pin2) annotation(
        Line(points = {{-70, -38}, {-104, -38}}, color = {0, 0, 255}));
      connect(resistor3.p, pin1) annotation(
        Line(points = {{-70, 12}, {-106, 12}, {-106, 10}}, color = {0, 0, 255}));
      connect(transformer.n1, idealTransformer.n1) annotation(
        Line(points = {{-12, 44}, {-80, 44}, {-80, -8}, {-10, -8}}, color = {0, 0, 255}));
      connect(idealTransformer1.n1, idealTransformer.n1) annotation(
        Line(points = {{-10, -58}, {-80, -58}, {-80, -8}, {-10, -8}}, color = {0, 0, 255}));
      connect(pin6, idealTransformer1.n1) annotation(
        Line(points = {{-104, -70}, {-10, -70}, {-10, -58}}, color = {0, 0, 255}));
    end Yn_d_Transformer;

    model COBEP_fault
      //filter parameters
      parameter Real Rf = 1.5649e-4;
      parameter Real Lf = 1.6604e-5;
      parameter Real Cf = 0.0013;
      parameter Real Rd = 0.03;
      // grid parameters
      parameter Real Rr = 3.27;
      parameter Real Lr = 1.23e-2;
      parameter Real Rcc = 0.1;
      parameter Real f = 60;
      //fault parameters
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
      Modelica.Electrical.Analog.Sensors.VoltageSensor vsa annotation(
        Placement(visible = true, transformation(origin = {213, 52}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.VoltageSensor vsb annotation(
        Placement(visible = true, transformation(origin = {209, -5}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.VoltageSensor vscv annotation(
        Placement(visible = true, transformation(origin = {214, -61}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.CurrentSensor ia annotation(
        Placement(visible = true, transformation(origin = {236, 33}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.CurrentSensor ib annotation(
        Placement(visible = true, transformation(origin = {238, -27}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.CurrentSensor ic annotation(
        Placement(visible = true, transformation(origin = {240, -87}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
      VSC_FZ.Testes.Clark3p clark3p1 annotation(
        Placement(visible = true, transformation(origin = {205, 146}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression(y = ia.i) annotation(
        Placement(visible = true, transformation(origin = {171, 144}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression1(y = ib.i) annotation(
        Placement(visible = true, transformation(origin = {171, 128}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression2(y = ic.i) annotation(
        Placement(visible = true, transformation(origin = {171, 112}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression4(y = vsb.v) annotation(
        Placement(visible = true, transformation(origin = {-129, 123}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression5(y = vscv.v) annotation(
        Placement(visible = true, transformation(origin = {-129, 107}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression6(y = vsa.v) annotation(
        Placement(visible = true, transformation(origin = {-129, 139}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Sources.SineVoltage sineVoltage(V = 11.394e3, freqHz = 60) annotation(
        Placement(visible = true, transformation(origin = {453, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Inductor inductor2(L = Lf) annotation(
        Placement(visible = true, transformation(origin = {120, -86}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Sources.SineVoltage sineVoltage1(V = 11.394e3, freqHz = 60, phase = -2.094395102393195) annotation(
        Placement(visible = true, transformation(origin = {455, -22}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Sources.SineVoltage sineVoltage2(V = 11.394e3, freqHz = 60, phase = 2.094395102393195) annotation(
        Placement(visible = true, transformation(origin = {455, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Ground ground1 annotation(
        Placement(visible = true, transformation(origin = {490, -22}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
      Modelica.Blocks.Sources.RealExpression realExpression7(y = dsogi_pll.d) annotation(
        Placement(visible = true, transformation(origin = {344, 173.5}, extent = {{-17, -9.5}, {17, 9.5}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression8(y = dsogi_pll.q) annotation(
        Placement(visible = true, transformation(origin = {311.5, 98.5}, extent = {{-15.5, -9.5}, {15.5, 9.5}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression9(y = Controle_potencia.idref) annotation(
        Placement(visible = true, transformation(origin = {333, 158.5}, extent = {{-28, -9.5}, {28, 9.5}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression10(y = Correntes.d) annotation(
        Placement(visible = true, transformation(origin = {341, 145}, extent = {{-19, -10}, {19, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression11(y = Correntes.q) annotation(
        Placement(visible = true, transformation(origin = {341.5, 126}, extent = {{-19.5, -10}, {19.5, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression12(y = Controle_potencia.iqref) annotation(
        Placement(visible = true, transformation(origin = {335, 113.5}, extent = {{-26, -8.5}, {26, 8.5}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression13(y = dsogi_pll.ang) annotation(
        Placement(visible = true, transformation(origin = {455, 101.5}, extent = {{-23, -11.5}, {23, 11.5}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant2(k = 0) annotation(
        Placement(visible = true, transformation(origin = {499, 102}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression14(y = clarkInv.a) annotation(
        Placement(visible = true, transformation(origin = {-73, -56.5}, extent = {{-17, -8.5}, {17, 8.5}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression15(y = clarkInv.b) annotation(
        Placement(visible = true, transformation(origin = {-75, -83.5}, extent = {{-17, -9.5}, {17, 9.5}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression16(y = clarkInv.c) annotation(
        Placement(visible = true, transformation(origin = {-76, -111.5}, extent = {{-18, -9.5}, {18, 9.5}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.VoltageSensor vdc_2 annotation(
        Placement(visible = true, transformation(origin = {-109, 21}, extent = {{-7, -7}, {7, 7}}, rotation = -90)));
      Modelica.Blocks.Math.Division division annotation(
        Placement(visible = true, transformation(origin = {444, 162}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Math.Division division1 annotation(
        Placement(visible = true, transformation(origin = {442, 125}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Nonlinear.Limiter limiter(limitsAtInit = true, uMax = 1, uMin = -1) annotation(
        Placement(visible = true, transformation(origin = {-28, -58}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
      Modelica.Blocks.Nonlinear.Limiter limiter1(limitsAtInit = true, uMax = 1, uMin = -1) annotation(
        Placement(visible = true, transformation(origin = {-28, -84}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
      Modelica.Blocks.Nonlinear.Limiter limiter2(limitsAtInit = true, uMax = 1, uMin = -1) annotation(
        Placement(visible = true, transformation(origin = {-28, -112}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
      Modelica.Blocks.Sources.Step Pref(height = -750e3, startTime = 0.4) annotation(
        Placement(visible = true, transformation(origin = {37, 140}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Step Qref(height = 0, startTime = 0) annotation(
        Placement(visible = true, transformation(origin = {37, 72}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      VSC_FZ.Testes.InversePark_fz inversePark_fz annotation(
        Placement(visible = true, transformation(origin = {484, 145}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      VSC_FZ.Testes.Park_Fz Correntes annotation(
        Placement(visible = true, transformation(origin = {243, 140}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Sources.ConstantVoltage constantVoltage1(V = 395.18) annotation(
        Placement(visible = true, transformation(origin = {-86, -16}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      Modelica.Electrical.Analog.Basic.Ground ground annotation(
        Placement(visible = true, transformation(origin = {-65, 70}, extent = {{-54, -86}, {-38, -70}}, rotation = 0)));
      Modelica.Blocks.Continuous.Filter filter(A_ripple = 0.1, analogFilter = Modelica.Blocks.Types.AnalogFilter.Bessel, f_cut = 60, filterType = Modelica.Blocks.Types.FilterType.LowPass) annotation(
        Placement(visible = true, transformation(origin = {497, -18}, extent = {{84, 4}, {104, 24}}, rotation = 0)));
      VSC_FZ.Testes.PowerCalc powerCalc annotation(
        Placement(visible = true, transformation(origin = {556, -11}, extent = {{-14, -14}, {14, 14}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression17(y = vsa.v) annotation(
        Placement(visible = true, transformation(origin = {522, 9}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression24(y = vsb.v) annotation(
        Placement(visible = true, transformation(origin = {522, -4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression25(y = vscv.v) annotation(
        Placement(visible = true, transformation(origin = {522, -17}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression26(y = ia.i) annotation(
        Placement(visible = true, transformation(origin = {522, -32}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression27(y = ib.i) annotation(
        Placement(visible = true, transformation(origin = {522, -48}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression28(y = ic.i) annotation(
        Placement(visible = true, transformation(origin = {522, -64}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Continuous.Filter filter1(A_ripple = 0.1, analogFilter = Modelica.Blocks.Types.AnalogFilter.Bessel, f_cut = 60, filterType = Modelica.Blocks.Types.FilterType.LowPass) annotation(
        Placement(visible = true, transformation(origin = {497, -50}, extent = {{84, 4}, {104, 24}}, rotation = 0)));
      VSC_FZ.Testes.PowerControl1 Controle_potencia(Pctrl_KI = 0.168, Pctrl_KP = 5.36e-4) annotation(
        Placement(visible = true, transformation(origin = {108.5, 104.5}, extent = {{-17.5, -17.5}, {17.5, 17.5}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression23(y = filter.y) annotation(
        Placement(visible = true, transformation(origin = {33.5, 111}, extent = {{-20.5, -10}, {20.5, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression29(y = filter1.y) annotation(
        Placement(visible = true, transformation(origin = {33.5, 94}, extent = {{-20.5, -10}, {20.5, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression18(y = vdc_2.v) annotation(
        Placement(visible = true, transformation(origin = {404, 78.5}, extent = {{-23, -11.5}, {23, 11.5}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Resistor resistor3(R = Rr, alpha = 0) annotation(
        Placement(visible = true, transformation(origin = {390, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Resistor resistor4(R = Rr, alpha = 0) annotation(
        Placement(visible = true, transformation(origin = {389, -22}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Resistor resistor5(R = Rr, alpha = 0) annotation(
        Placement(visible = true, transformation(origin = {389, -82}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Ground ground2 annotation(
        Placement(visible = true, transformation(origin = {229.5, -5.5}, extent = {{-6.5, -6.5}, {6.5, 6.5}}, rotation = 90)));
      Modelica.Electrical.Analog.Basic.Ground ground3 annotation(
        Placement(visible = true, transformation(origin = {231.5, 51.5}, extent = {{-6.5, -6.5}, {6.5, 6.5}}, rotation = 90)));
      Modelica.Electrical.Analog.Basic.Ground ground4 annotation(
        Placement(visible = true, transformation(origin = {236.5, -61.5}, extent = {{-6.5, -6.5}, {6.5, 6.5}}, rotation = 90)));
      Modelica.Electrical.Analog.Basic.Inductor inductor3(L = Lr) annotation(
        Placement(visible = true, transformation(origin = {422, 39}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Inductor inductor4(L = Lr) annotation(
        Placement(visible = true, transformation(origin = {421, -22}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Inductor inductor5(L = Lr) annotation(
        Placement(visible = true, transformation(origin = {421, -82}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Ground ground5 annotation(
        Placement(visible = true, transformation(origin = {393.5, 13.5}, extent = {{-6.5, -6.5}, {6.5, 6.5}}, rotation = 90)));
      Modelica.Electrical.Analog.Sensors.VoltageSensor pca annotation(
        Placement(visible = true, transformation(origin = {375, 14}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.VoltageSensor pcc annotation(
        Placement(visible = true, transformation(origin = {373, -107}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Ground ground6 annotation(
        Placement(visible = true, transformation(origin = {391.5, -107.5}, extent = {{-6.5, -6.5}, {6.5, 6.5}}, rotation = 90)));
      Modelica.Electrical.Analog.Basic.Ground ground7 annotation(
        Placement(visible = true, transformation(origin = {396.5, -45.5}, extent = {{-6.5, -6.5}, {6.5, 6.5}}, rotation = 90)));
      Modelica.Electrical.Analog.Sensors.VoltageSensor pcb annotation(
        Placement(visible = true, transformation(origin = {378, -45}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
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
        Placement(visible = true, transformation(origin = {597, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Division P_pu annotation(
        Placement(visible = true, transformation(origin = {630, 46}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant7(k = 915e3) annotation(
        Placement(visible = true, transformation(origin = {597, -89}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Division Q_pu annotation(
        Placement(visible = true, transformation(origin = {630, -73}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant8(k = 1312.159) annotation(
        Placement(visible = true, transformation(origin = {262, 175}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
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
      Modelica.Electrical.Analog.Sensors.CurrentSensor ipca annotation(
        Placement(visible = true, transformation(origin = {351, 38}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.CurrentSensor ipcb annotation(
        Placement(visible = true, transformation(origin = {353, -22}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.CurrentSensor ipcc annotation(
        Placement(visible = true, transformation(origin = {355, -82}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression39(y = ia.i) annotation(
        Placement(visible = true, transformation(origin = {-126, 219}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression40(y = ib.i) annotation(
        Placement(visible = true, transformation(origin = {-126, 203}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      VSC_FZ.Testes.PLL_fz2 correntes_tot(frequency = 60) annotation(
        Placement(visible = true, transformation(origin = {-26, 217}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant14(k = 1698) annotation(
        Placement(visible = true, transformation(origin = {-47, 164}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Division itot_d annotation(
        Placement(visible = true, transformation(origin = {-7, 189}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Math.Division itot_q annotation(
        Placement(visible = true, transformation(origin = {-8, 165}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      VSC_FZ.Testes.Clark3p clark3p2 annotation(
        Placement(visible = true, transformation(origin = {-79, 213}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression38(y = ic.i) annotation(
        Placement(visible = true, transformation(origin = {-126, 187}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant1(k = 915e3) annotation(
        Placement(visible = true, transformation(origin = {746, 33}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression22(y = ipcb.i) annotation(
        Placement(visible = true, transformation(origin = {671, -45}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant13(k = 915e3) annotation(
        Placement(visible = true, transformation(origin = {746, -86}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Continuous.Filter filter2(analogFilter = Modelica.Blocks.Types.AnalogFilter.Bessel, f_cut = 60, filterType = Modelica.Blocks.Types.FilterType.LowPass) annotation(
        Placement(visible = true, transformation(origin = {646, -47}, extent = {{84, 4}, {104, 24}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression33(y = ipca.i) annotation(
        Placement(visible = true, transformation(origin = {671, -29}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Division P_pupc annotation(
        Placement(visible = true, transformation(origin = {779, 49}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Continuous.Filter filter3(analogFilter = Modelica.Blocks.Types.AnalogFilter.Bessel, f_cut = 60, filterType = Modelica.Blocks.Types.FilterType.LowPass) annotation(
        Placement(visible = true, transformation(origin = {646, -16}, extent = {{84, 4}, {104, 24}}, rotation = 0)));
      Modelica.Blocks.Math.Division Q_pupc annotation(
        Placement(visible = true, transformation(origin = {779, -70}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression34(y = pcb.v) annotation(
        Placement(visible = true, transformation(origin = {671, -1}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression35(y = ipcc.i) annotation(
        Placement(visible = true, transformation(origin = {671, -61}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression36(y = pcc.v) annotation(
        Placement(visible = true, transformation(origin = {671, -14}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      VSC_FZ.Testes.PowerCalc powerCalc_pcc annotation(
        Placement(visible = true, transformation(origin = {705, -9}, extent = {{-14, -14}, {14, 14}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression37(y = pca.v) annotation(
        Placement(visible = true, transformation(origin = {671, 12}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
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
      Modelica.Blocks.Sources.Constant constant15(k = 41.837) annotation(
        Placement(visible = true, transformation(origin = {630, 217}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression43(y = ipca.i) annotation(
        Placement(visible = true, transformation(origin = {511, 279}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Division ipq annotation(
        Placement(visible = true, transformation(origin = {669, 218}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Electrical.Analog.Ideal.IdealClosingSwitch idealClosingSwitch annotation(
        Placement(visible = true, transformation(origin = {428, -48}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
      Modelica.Electrical.Analog.Basic.Resistor resistor6(R = Rcc, alpha = 0) annotation(
        Placement(visible = true, transformation(origin = {443, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Ground ground11 annotation(
        Placement(visible = true, transformation(origin = {463.5, 0.5}, extent = {{-6.5, -6.5}, {6.5, 6.5}}, rotation = 90)));
      Modelica.Electrical.Analog.Ideal.IdealClosingSwitch idealClosingSwitch1 annotation(
        Placement(visible = true, transformation(origin = {428, -106}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
      Modelica.Electrical.Analog.Ideal.IdealClosingSwitch switch2 annotation(
        Placement(visible = true, transformation(origin = {425, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
      Modelica.Blocks.MathBoolean.And and1(nu = 2) annotation(
        Placement(visible = true, transformation(origin = {385, -156}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
      Modelica.Blocks.Sources.BooleanStep booleanStep4(startTime = 0, startValue = true) annotation(
        Placement(visible = true, transformation(origin = {358, -141}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Ground ground13 annotation(
        Placement(visible = true, transformation(origin = {465.5, -115.5}, extent = {{-6.5, -6.5}, {6.5, 6.5}}, rotation = 90)));
      Modelica.Electrical.Analog.Basic.Resistor resistor7(R = Rcc, alpha = 0) annotation(
        Placement(visible = true, transformation(origin = {445, -58}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Ground ground12 annotation(
        Placement(visible = true, transformation(origin = {465.5, -57.5}, extent = {{-6.5, -6.5}, {6.5, 6.5}}, rotation = 90)));
      Modelica.Blocks.Sources.BooleanStep booleanStep3(startTime = 0, startValue = false) annotation(
        Placement(visible = true, transformation(origin = {357, -174}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Resistor resistor8(R = Rcc, alpha = 0) annotation(
        Placement(visible = true, transformation(origin = {445, -118}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.BooleanStep booleanStep(startTime = 1, startValue = false) annotation(
        Placement(visible = true, transformation(origin = {285, -141}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.MathBoolean.And and2(nu = 2) annotation(
        Placement(visible = true, transformation(origin = {311, -156}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
      Modelica.Blocks.Sources.BooleanStep booleanStep1(startTime = 0, startValue = false) annotation(
        Placement(visible = true, transformation(origin = {283, -174}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.BooleanStep booleanStep2(startTime = 0, startValue = true) annotation(
        Placement(visible = true, transformation(origin = {404, -186}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.MathBoolean.And and3(nu = 2) annotation(
        Placement(visible = true, transformation(origin = {430, -201}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
      Modelica.Blocks.Sources.BooleanStep booleanStep5(startTime = 0, startValue = false) annotation(
        Placement(visible = true, transformation(origin = {402, -219}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
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
      VSC_FZ.Testes.ClarkInv clarkInv annotation(
        Placement(visible = true, transformation(origin = {550, 149}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Capacitor capacitor(C = Cf) annotation(
        Placement(visible = true, transformation(origin = {133, 20}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      Modelica.Electrical.Analog.Basic.Resistor resistor9(R = Rd, alpha = 0) annotation(
        Placement(visible = true, transformation(origin = {154, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Ground ground14 annotation(
        Placement(visible = true, transformation(origin = {175.5, 9.5}, extent = {{-6.5, -6.5}, {6.5, 6.5}}, rotation = 90)));
      Modelica.Electrical.Analog.Basic.Resistor resistor10(R = Rf, alpha = 0) annotation(
        Placement(visible = true, transformation(origin = {158, 34}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Inductor inductor6(L = Lf) annotation(
        Placement(visible = true, transformation(origin = {186, 34}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Ground ground15 annotation(
        Placement(visible = true, transformation(origin = {176.5, -50.5}, extent = {{-6.5, -6.5}, {6.5, 6.5}}, rotation = 90)));
      Modelica.Electrical.Analog.Basic.Resistor resistor11(R = Rf, alpha = 0) annotation(
        Placement(visible = true, transformation(origin = {159, -25}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Capacitor capacitor1(C = Cf) annotation(
        Placement(visible = true, transformation(origin = {135, -40}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      Modelica.Electrical.Analog.Basic.Inductor inductor7(L = Lf) annotation(
        Placement(visible = true, transformation(origin = {187, -25}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Resistor resistor12(R = Rd, alpha = 0) annotation(
        Placement(visible = true, transformation(origin = {155, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Ground ground16 annotation(
        Placement(visible = true, transformation(origin = {179.5, -110.5}, extent = {{-6.5, -6.5}, {6.5, 6.5}}, rotation = 90)));
      Modelica.Electrical.Analog.Basic.Resistor resistor13(R = Rf, alpha = 0) annotation(
        Placement(visible = true, transformation(origin = {162, -86}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Capacitor capacitor2(C = Cf) annotation(
        Placement(visible = true, transformation(origin = {137, -100}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      Modelica.Electrical.Analog.Basic.Inductor inductor8(L = Lf) annotation(
        Placement(visible = true, transformation(origin = {190, -86}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Resistor resistor14(R = Rd, alpha = 0) annotation(
        Placement(visible = true, transformation(origin = {158, -110}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      VSC_FZ.Testes.DSOGI_PLL dsogi_pll annotation(
        Placement(visible = true, transformation(origin = {-81, 124}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      VSC_FZ.Testes.InnerControl_PIR innerControl_PIR(Ictrl_KI = 0.4916, Ictrl_KP = 5.216e-2, KIh = 600, freq = 2*3.1416*60, h = 2, lw = 2*3.1416*60*2*Lf, r = 2*Rf) annotation(
        Placement(visible = true, transformation(origin = {402.5, 150.5}, extent = {{-15.5, -15.5}, {15.5, 15.5}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression3(y = dsogi_pll.ang) annotation(
        Placement(visible = true, transformation(origin = {211.5, 118}, extent = {{-20.5, -9}, {20.5, 9}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression30(y = pcc.v) annotation(
        Placement(visible = true, transformation(origin = {77, 190}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression31(y = pcb.v) annotation(
        Placement(visible = true, transformation(origin = {77, 206}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression32(y = pca.v) annotation(
        Placement(visible = true, transformation(origin = {77, 222}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Clark3p clark3p3 annotation(
        Placement(visible = true, transformation(origin = {121, 209}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Division pcq annotation(
        Placement(visible = true, transformation(origin = {212, 186}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant3(k = 11.26765e3) annotation(
        Placement(visible = true, transformation(origin = {173, 186}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      PLL_fz2 pLL_fz21(frequency = 60) annotation(
        Placement(visible = true, transformation(origin = {174, 213}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Division pcd annotation(
        Placement(visible = true, transformation(origin = {213, 210}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      VSC_FZ.Testes.D_D_Transformer d_D_Transformer(L_trans_prim = 3.06e-5, L_trans_secon = 3.01e-3, R_trans_prim = 11.56e-4, R_trans_secon = 1.137, Vac_primary = 440, Vac_secondary = 13800) annotation(
        Placement(visible = true, transformation(origin = {299, -24}, extent = {{-27, -27}, {27, 27}}, rotation = 0)));
      Modelica.Blocks.Math.Division division2 annotation(
        Placement(visible = true, transformation(origin = {58, 313}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Math.Division division3 annotation(
        Placement(visible = true, transformation(origin = {58, 285}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression47(y = ipcb.i) annotation(
        Placement(visible = true, transformation(origin = {1, 287}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Division division4 annotation(
        Placement(visible = true, transformation(origin = {58, 259}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression48(y = ipcc.i) annotation(
        Placement(visible = true, transformation(origin = {1, 271}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression49(y = ipca.i) annotation(
        Placement(visible = true, transformation(origin = {1, 303}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant17(k = 59.166) annotation(
        Placement(visible = true, transformation(origin = {21, 240}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Ramp ramp(duration = 0.331, height = 750e3, startTime = 0.4) annotation(
        Placement(visible = true, transformation(origin = {38, 171}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Ramp ramp1(duration = 0, height = 0, startTime = 0) annotation(
        Placement(visible = true, transformation(origin = {75, 68}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant18(k = 915e3) annotation(
        Placement(visible = true, transformation(origin = {114, 61}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Division Qref_pu annotation(
        Placement(visible = true, transformation(origin = {147, 76}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Math.Division pca_pu annotation(
        Placement(visible = true, transformation(origin = {627, 112}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant12(k = 11.26765e3) annotation(
        Placement(visible = true, transformation(origin = {594, 99}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression20(y = pca.v) annotation(
        Placement(visible = true, transformation(origin = {597, 119}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(constantVoltage.p, vsc.pin_p) annotation(
        Line(points = {{-86, 28}, {-86, 34}, {-15, 34}, {-15, 33.5}}, color = {0, 0, 255}));
      connect(resistor.n, inductor.p) annotation(
        Line(points = {{97, 34}, {105, 34}}, color = {0, 0, 255}));
      connect(resistor1.n, inductor1.p) annotation(
        Line(points = {{100, -26}, {108, -26}}, color = {0, 0, 255}));
      connect(realExpression.y, clark3p1.V_abc_D[1]) annotation(
        Line(points = {{182, 144}, {187, 144}, {187, 146}, {193, 146}}, color = {0, 0, 127}));
      connect(realExpression1.y, clark3p1.V_abc_D[2]) annotation(
        Line(points = {{182, 128}, {184, 128}, {184, 146}, {193, 146}}, color = {0, 0, 127}));
      connect(realExpression2.y, clark3p1.V_abc_D[3]) annotation(
        Line(points = {{182, 112}, {184, 112}, {184, 146}, {193, 146}}, color = {0, 0, 127}));
      connect(inductor2.p, resistor2.n) annotation(
        Line(points = {{110, -86}, {102, -86}}, color = {0, 0, 255}));
      connect(sineVoltage.n, ground1.p) annotation(
        Line(points = {{463, 40}, {463, 38}, {480, 38}, {480, -22}}, color = {0, 0, 255}));
      connect(sineVoltage1.n, ground1.p) annotation(
        Line(points = {{465, -22}, {480, -22}}, color = {0, 0, 255}));
      connect(sineVoltage2.n, ground1.p) annotation(
        Line(points = {{465, -80}, {465, -82}, {480, -82}, {480, -22}}, color = {0, 0, 255}));
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
      connect(division.y, inversePark_fz.d) annotation(
        Line(points = {{452, 162}, {463.7, 162}, {463.7, 149}, {471.7, 149}}, color = {0, 0, 127}));
      connect(division1.y, inversePark_fz.q) annotation(
        Line(points = {{449.7, 125}, {465.7, 125}, {465.7, 141}, {471.7, 141}}, color = {0, 0, 127}));
      connect(inversePark_fz.theta, realExpression13.y) annotation(
        Line(points = {{484, 133}, {484, 101.5}, {480, 101.5}}, color = {0, 0, 127}));
      connect(clark3p1.V_ab_D[1], Correntes.alpha) annotation(
        Line(points = {{216, 146}, {223.5, 146}, {223.5, 144}, {231, 144}}, color = {0, 0, 127}));
      connect(clark3p1.V_ab_D[2], Correntes.beta) annotation(
        Line(points = {{216, 146}, {224, 146}, {224, 136}, {231, 136}}, color = {0, 0, 127}));
      connect(constantVoltage1.n, vsc.pin_n) annotation(
        Line(points = {{-86, -26}, {-60, -26}, {-60, 20.5}, {-15, 20.5}}, color = {0, 0, 255}));
      connect(constantVoltage.n, constantVoltage1.p) annotation(
        Line(points = {{-86, 8}, {-86, -6}}, color = {0, 0, 255}));
      connect(realExpression17.y, powerCalc.V1) annotation(
        Line(points = {{533, 9}, {538, 9}, {538, 1}, {541, 1}}, color = {0, 0, 127}));
      connect(realExpression24.y, powerCalc.V2) annotation(
        Line(points = {{533, -4}, {537, -4}, {537, -3}, {541, -3}}, color = {0, 0, 127}));
      connect(realExpression25.y, powerCalc.V3) annotation(
        Line(points = {{533, -17}, {535, -17}, {535, -7}, {541, -7}}, color = {0, 0, 127}));
      connect(realExpression26.y, powerCalc.I1) annotation(
        Line(points = {{533, -32}, {537, -32}, {537, -14}, {541, -14}}, color = {0, 0, 127}));
      connect(realExpression27.y, powerCalc.I2) annotation(
        Line(points = {{533, -48}, {538, -48}, {538, -18}, {541, -18}}, color = {0, 0, 127}));
      connect(realExpression28.y, powerCalc.I3) annotation(
        Line(points = {{533, -64}, {539, -64}, {539, -23}, {541, -23}}, color = {0, 0, 127}));
      connect(powerCalc.P, filter.u) annotation(
        Line(points = {{571, -5}, {571, -4}, {579, -4}}, color = {0, 0, 127}));
      connect(powerCalc.Q, filter1.u) annotation(
        Line(points = {{571, -15}, {574.4, -15}, {574.4, -36.2}, {579.4, -36.2}}, color = {0, 0, 127}));
      connect(realExpression29.y, Controle_potencia.Q) annotation(
        Line(points = {{56, 94}, {72, 94}, {72, 100}, {89, 100}}, color = {0, 0, 127}));
      connect(realExpression23.y, Controle_potencia.P) annotation(
        Line(points = {{56, 111}, {67, 111}, {67, 106}, {89, 106}}, color = {0, 0, 127}));
      connect(realExpression18.y, division1.u2) annotation(
        Line(points = {{429, 78.5}, {429, 121}, {434, 121}}, color = {0, 0, 127}));
      connect(realExpression18.y, division.u2) annotation(
        Line(points = {{429, 78.5}, {429, 158}, {436, 158}}, color = {0, 0, 127}));
      connect(vsb.n, ground2.p) annotation(
        Line(points = {{216, -5}, {223, -5}}, color = {0, 0, 255}));
      connect(vsa.n, ground3.p) annotation(
        Line(points = {{220, 52}, {225, 52}}, color = {0, 0, 255}));
      connect(vscv.n, ground4.p) annotation(
        Line(points = {{222, -61}, {230, -61}}, color = {0, 0, 255}));
      connect(resistor4.n, inductor4.n) annotation(
        Line(points = {{399, -22}, {411, -22}}, color = {0, 0, 255}));
      connect(resistor5.n, inductor5.n) annotation(
        Line(points = {{399, -82}, {411, -82}}, color = {0, 0, 255}));
      connect(inductor5.p, sineVoltage2.p) annotation(
        Line(points = {{431, -82}, {438, -82}, {438, -80}, {445, -80}}, color = {0, 0, 255}));
      connect(inductor4.p, sineVoltage1.p) annotation(
        Line(points = {{431, -22}, {445, -22}}, color = {0, 0, 255}));
      connect(resistor3.n, inductor3.p) annotation(
        Line(points = {{400, 38}, {404, 38}, {404, 39}, {412, 39}}, color = {0, 0, 255}));
      connect(inductor3.n, sineVoltage.p) annotation(
        Line(points = {{432, 39}, {432, 40}, {443, 40}}, color = {0, 0, 255}));
      connect(pca.n, ground5.p) annotation(
        Line(points = {{382, 14}, {387, 14}}, color = {0, 0, 255}));
      connect(pcc.n, ground6.p) annotation(
        Line(points = {{380, -107}, {385, -107}}, color = {0, 0, 255}));
      connect(pcb.n, ground7.p) annotation(
        Line(points = {{385, -45}, {390, -45}}, color = {0, 0, 255}));
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
      connect(constant4.y, vsq.u2) annotation(
        Line(points = {{-44, 80}, {-24, 80}, {-24, 77}}, color = {0, 0, 127}));
      connect(constant4.y, vsd.u2) annotation(
        Line(points = {{-44, 80}, {-31, 80}, {-31, 101}, {-23, 101}}, color = {0, 0, 127}));
      connect(constant5.y, Pref_pu.u2) annotation(
        Line(points = {{108, 143}, {113, 143}, {113, 155}, {122, 155}}, color = {0, 0, 127}));
      connect(constant6.y, P_pu.u2) annotation(
        Line(points = {{608, 30}, {613, 30}, {613, 42}, {622, 42}}, color = {0, 0, 127}));
      connect(powerCalc.P, P_pu.u1) annotation(
        Line(points = {{571, -5}, {573.4, -5}, {573.4, 50.44}, {622.4, 50.44}}, color = {0, 0, 127}));
      connect(constant7.y, Q_pu.u2) annotation(
        Line(points = {{608, -89}, {613, -89}, {613, -77}, {622, -77}}, color = {0, 0, 127}));
      connect(powerCalc.Q, Q_pu.u1) annotation(
        Line(points = {{571, -15}, {575.4, -15}, {575.4, -69.2}, {622.4, -69.2}}, color = {0, 0, 127}));
      connect(constant8.y, id_pu.u2) annotation(
        Line(points = {{273, 175}, {273, 185}, {288, 185}}, color = {0, 0, 127}));
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
      connect(resistor3.p, ipca.n) annotation(
        Line(points = {{380, 38}, {357, 38}}, color = {0, 0, 255}));
      connect(resistor4.p, ipcb.n) annotation(
        Line(points = {{379, -22}, {359, -22}}, color = {0, 0, 255}));
      connect(resistor5.p, ipcc.n) annotation(
        Line(points = {{379, -82}, {361, -82}}, color = {0, 0, 255}));
      connect(constant14.y, itot_q.u2) annotation(
        Line(points = {{-36, 164}, {-16, 164}, {-16, 161}}, color = {0, 0, 127}));
      connect(realExpression39.y, clark3p2.V_abc_D[1]) annotation(
        Line(points = {{-115, 219}, {-99, 219}, {-99, 213}, {-91, 213}}, color = {0, 0, 127}));
      connect(clark3p2.V_ab_D[1], correntes_tot.alpha) annotation(
        Line(points = {{-68, 213}, {-55.5, 213}, {-55.5, 221}, {-38, 221}}, color = {0, 0, 127}));
      connect(correntes_tot.q, itot_q.u1) annotation(
        Line(points = {{-15, 208}, {-29, 208}, {-29, 166.4}, {-16, 166.4}}, color = {0, 0, 127}));
      connect(constant14.y, itot_d.u2) annotation(
        Line(points = {{-36, 164}, {-23, 164}, {-23, 185}, {-15, 185}}, color = {0, 0, 127}));
      connect(realExpression40.y, clark3p2.V_abc_D[2]) annotation(
        Line(points = {{-115, 203}, {-99, 203}, {-99, 213}, {-91, 213}}, color = {0, 0, 127}));
      connect(correntes_tot.beta, clark3p2.V_ab_D[2]) annotation(
        Line(points = {{-38.2, 212.4}, {-51.2, 212.4}, {-51.2, 213.4}, {-68.2, 213.4}}, color = {0, 0, 127}));
      connect(realExpression38.y, clark3p2.V_abc_D[3]) annotation(
        Line(points = {{-115, 187}, {-99, 187}, {-99, 213}, {-91, 213}}, color = {0, 0, 127}));
      connect(correntes_tot.d, itot_d.u1) annotation(
        Line(points = {{-15, 214}, {-25, 214}, {-25, 190}, {-15, 190}}, color = {0, 0, 127}));
      connect(constant1.y, P_pupc.u2) annotation(
        Line(points = {{757, 33}, {762, 33}, {762, 45}, {771, 45}}, color = {0, 0, 127}));
      connect(realExpression33.y, powerCalc_pcc.I1) annotation(
        Line(points = {{682, -29}, {686, -29}, {686, -12}, {690, -12}}, color = {0, 0, 127}));
      connect(realExpression35.y, powerCalc_pcc.I3) annotation(
        Line(points = {{682, -61}, {688, -61}, {688, -21}, {690, -21}}, color = {0, 0, 127}));
      connect(constant13.y, Q_pupc.u2) annotation(
        Line(points = {{757, -86}, {762, -86}, {762, -74}, {771, -74}}, color = {0, 0, 127}));
      connect(realExpression37.y, powerCalc_pcc.V1) annotation(
        Line(points = {{682, 12}, {687, 12}, {687, 3}, {690, 3}}, color = {0, 0, 127}));
      connect(powerCalc_pcc.Q, Q_pupc.u1) annotation(
        Line(points = {{720.4, -13.2}, {724.4, -13.2}, {724.4, -66.2}, {771.4, -66.2}}, color = {0, 0, 127}));
      connect(realExpression22.y, powerCalc_pcc.I2) annotation(
        Line(points = {{682, -45}, {687, -45}, {687, -16}, {690, -16}}, color = {0, 0, 127}));
      connect(realExpression34.y, powerCalc_pcc.V2) annotation(
        Line(points = {{682, -1}, {690, -1}}, color = {0, 0, 127}));
      connect(realExpression36.y, powerCalc_pcc.V3) annotation(
        Line(points = {{682, -14}, {684, -14}, {684, -5}, {690, -5}}, color = {0, 0, 127}));
      connect(powerCalc_pcc.Q, filter2.u) annotation(
        Line(points = {{720.4, -13.2}, {723.4, -13.2}, {723.4, -33.2}, {728.4, -33.2}}, color = {0, 0, 127}));
      connect(powerCalc_pcc.P, filter3.u) annotation(
        Line(points = {{720.4, -2.56}, {728.4, -2.56}, {728.4, -1.56}}, color = {0, 0, 127}));
      connect(powerCalc_pcc.P, P_pupc.u1) annotation(
        Line(points = {{720.4, -2.56}, {722.4, -2.56}, {722.4, 53.44}, {771.4, 53.44}}, color = {0, 0, 127}));
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
        Line(points = {{368, -174}, {376, -174}, {376, -156}, {379, -156}}, color = {255, 0, 255}));
      connect(booleanStep4.y, and1.u[1]) annotation(
        Line(points = {{369, -141}, {376, -141}, {376, -156}, {379, -156}}, color = {255, 0, 255}));
      connect(switch2.n, resistor3.p) annotation(
        Line(points = {{425, 20}, {425, 25}, {380, 25}, {380, 38}}, color = {0, 0, 255}));
      connect(resistor6.p, switch2.p) annotation(
        Line(points = {{433, 0}, {425, 0}}, color = {0, 0, 255}));
      connect(ground11.p, resistor6.n) annotation(
        Line(points = {{457, 0.5}, {453, 0.5}, {453, -0.5}}, color = {0, 0, 255}));
      connect(idealClosingSwitch.n, resistor4.p) annotation(
        Line(points = {{428, -38}, {379, -38}, {379, -22}}, color = {0, 0, 255}));
      connect(idealClosingSwitch.p, resistor7.p) annotation(
        Line(points = {{428, -58}, {435, -58}}, color = {0, 0, 255}));
      connect(resistor7.n, ground12.p) annotation(
        Line(points = {{455, -58}, {459, -58}, {459, -57}}, color = {0, 0, 255}));
      connect(idealClosingSwitch1.n, resistor5.p) annotation(
        Line(points = {{428, -96}, {375, -96}, {375, -82}, {379, -82}}, color = {0, 0, 255}));
      connect(idealClosingSwitch1.p, resistor8.p) annotation(
        Line(points = {{428, -116}, {435, -116}, {435, -118}}, color = {0, 0, 255}));
      connect(resistor8.n, ground13.p) annotation(
        Line(points = {{455, -118}, {459, -118}, {459, -115}}, color = {0, 0, 255}));
      connect(booleanStep.y, and2.u[1]) annotation(
        Line(points = {{296, -141}, {302, -141}, {302, -156}, {305, -156}}, color = {255, 0, 255}));
      connect(booleanStep1.y, and2.u[2]) annotation(
        Line(points = {{294, -174}, {302, -174}, {302, -156}, {305, -156}}, color = {255, 0, 255}));
      connect(booleanStep2.y, and3.u[1]) annotation(
        Line(points = {{415, -186}, {421, -186}, {421, -201}, {424, -201}}, color = {255, 0, 255}));
      connect(booleanStep5.y, and3.u[2]) annotation(
        Line(points = {{413, -219}, {421, -219}, {421, -201}, {424, -201}}, color = {255, 0, 255}));
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
        Line(points = {{495, 149}, {516, 149}, {516, 154}, {537, 154}}, color = {0, 0, 127}));
      connect(inversePark_fz.beta, clarkInv.B) annotation(
        Line(points = {{495, 141}, {524, 141}, {524, 148}, {537, 148}}, color = {0, 0, 127}));
      connect(constant2.y, clarkInv.C) annotation(
        Line(points = {{510, 102}, {530, 102}, {530, 143}, {537, 143}}, color = {0, 0, 127}));
      connect(capacitor.n, resistor9.p) annotation(
        Line(points = {{133, 10}, {144, 10}}, color = {0, 0, 255}));
      connect(capacitor.p, inductor.n) annotation(
        Line(points = {{133, 30}, {133, 34}, {125, 34}}, color = {0, 0, 255}));
      connect(resistor9.n, ground14.p) annotation(
        Line(points = {{164, 10}, {166.5, 10}, {166.5, 9.5}, {169, 9.5}}, color = {0, 0, 255}));
      connect(capacitor.p, resistor10.p) annotation(
        Line(points = {{133, 30}, {133, 34}, {148, 34}}, color = {0, 0, 255}));
      connect(resistor10.n, inductor6.p) annotation(
        Line(points = {{168, 34}, {176, 34}}, color = {0, 0, 255}));
      connect(inductor1.n, resistor11.p) annotation(
        Line(points = {{128, -26}, {149, -26}, {149, -25}}, color = {0, 0, 255}));
      connect(capacitor1.p, inductor1.n) annotation(
        Line(points = {{135, -30}, {135, -26}, {128, -26}}, color = {0, 0, 255}));
      connect(capacitor1.n, resistor12.p) annotation(
        Line(points = {{135, -50}, {145, -50}}, color = {0, 0, 255}));
      connect(resistor11.n, inductor7.p) annotation(
        Line(points = {{169, -25}, {177, -25}}, color = {0, 0, 255}));
      connect(resistor12.n, ground15.p) annotation(
        Line(points = {{165, -50}, {170, -50}}, color = {0, 0, 255}));
      connect(inductor2.n, resistor13.p) annotation(
        Line(points = {{130, -86}, {152, -86}}, color = {0, 0, 255}));
      connect(capacitor2.p, inductor2.n) annotation(
        Line(points = {{137, -90}, {137, -86}, {130, -86}}, color = {0, 0, 255}));
      connect(capacitor2.n, resistor14.p) annotation(
        Line(points = {{137, -110}, {148, -110}}, color = {0, 0, 255}));
      connect(resistor13.n, inductor8.p) annotation(
        Line(points = {{172, -86}, {180, -86}}, color = {0, 0, 255}));
      connect(resistor14.n, ground16.p) annotation(
        Line(points = {{168, -110}, {173, -110}}, color = {0, 0, 255}));
      connect(and2.y, switch2.control) annotation(
        Line(points = {{318, -156}, {324, -156}, {324, -121}, {399, -121}, {399, 10}, {413, 10}}, color = {255, 0, 255}));
      connect(and1.y, idealClosingSwitch.control) annotation(
        Line(points = {{392, -156}, {405, -156}, {405, -48}, {416, -48}}, color = {255, 0, 255}));
      connect(and3.y, idealClosingSwitch1.control) annotation(
        Line(points = {{437, -201}, {449, -201}, {449, -141}, {412, -141}, {412, -106}, {416, -106}}, color = {255, 0, 255}));
      connect(realExpression6.y, dsogi_pll.V_abc[1]) annotation(
        Line(points = {{-118, 139}, {-109, 139}, {-109, 124}, {-92, 124}}, color = {0, 0, 127}));
      connect(realExpression4.y, dsogi_pll.V_abc[2]) annotation(
        Line(points = {{-118, 123}, {-105, 123}, {-105, 124}, {-92, 124}}, color = {0, 0, 127}));
      connect(realExpression5.y, dsogi_pll.V_abc[3]) annotation(
        Line(points = {{-118, 107}, {-109, 107}, {-109, 124}, {-92, 124}}, color = {0, 0, 127}));
      connect(dsogi_pll.d, vsd.u1) annotation(
        Line(points = {{-70, 133}, {-56, 133}, {-56, 109}, {-23, 109}}, color = {0, 0, 127}));
      connect(dsogi_pll.q, vsq.u1) annotation(
        Line(points = {{-70, 127}, {-37, 127}, {-37, 85}, {-24, 85}}, color = {0, 0, 127}));
      connect(realExpression7.y, innerControl_PIR.vd) annotation(
        Line(points = {{363, 174}, {363, 175}, {385, 175}, {385, 164}}, color = {0, 0, 127}));
      connect(realExpression10.y, innerControl_PIR.Id) annotation(
        Line(points = {{362, 145}, {370.5, 145}, {370.5, 152}, {385, 152}}, color = {0, 0, 127}));
      connect(realExpression11.y, innerControl_PIR.Iq) annotation(
        Line(points = {{363, 126}, {373.5, 126}, {373.5, 146}, {385, 146}}, color = {0, 0, 127}));
      connect(innerControl_PIR.vdref, division.u1) annotation(
        Line(points = {{420, 159}, {421.5, 159}, {421.5, 166}, {436, 166}}, color = {0, 0, 127}));
      connect(innerControl_PIR.vqref, division1.u1) annotation(
        Line(points = {{420, 143}, {420, 129}, {434, 129}}, color = {0, 0, 127}));
      connect(innerControl_PIR.vq, realExpression8.y) annotation(
        Line(points = {{385, 137}, {380, 137}, {380, 99}, {329, 99}}, color = {0, 0, 127}));
      connect(realExpression3.y, Correntes.theta) annotation(
        Line(points = {{234, 118}, {243, 118}, {243, 128}}, color = {0, 0, 127}));
      connect(vsa.p, inductor6.n) annotation(
        Line(points = {{206, 52}, {196, 52}, {196, 34}}, color = {0, 0, 255}));
      connect(vsb.p, inductor7.n) annotation(
        Line(points = {{202, -5}, {197, -5}, {197, -25}}, color = {0, 0, 255}));
      connect(vscv.p, inductor8.n) annotation(
        Line(points = {{206, -61}, {200, -61}, {200, -86}}, color = {0, 0, 255}));
      connect(vsc.Vta, resistor.p) annotation(
        Line(points = {{13, 35}, {45, 35}, {45, 34}, {77, 34}}, color = {0, 0, 255}));
      connect(vsc.Vtb, resistor1.p) annotation(
        Line(points = {{13, 27}, {54, 27}, {54, -26}, {80, -26}}, color = {0, 0, 255}));
      connect(vsc.Vtc, resistor2.p) annotation(
        Line(points = {{13, 19}, {44, 19}, {44, -86}, {82, -86}}, color = {0, 0, 255}));
      connect(inductor6.n, ia.p) annotation(
        Line(points = {{196, 34}, {230, 34}, {230, 33}}, color = {0, 0, 255}));
      connect(inductor7.n, ib.p) annotation(
        Line(points = {{197, -25}, {232, -25}, {232, -27}}, color = {0, 0, 255}));
      connect(inductor8.n, ic.p) annotation(
        Line(points = {{200, -86}, {216.5, -86}, {216.5, -87}, {234, -87}}, color = {0, 0, 255}));
      connect(realExpression9.y, innerControl_PIR.Id_ref) annotation(
        Line(points = {{364, 159}, {385, 159}, {385, 158}}, color = {0, 0, 127}));
      connect(realExpression12.y, innerControl_PIR.Iqref) annotation(
        Line(points = {{364, 114}, {377, 114}, {377, 141}, {385, 141}}, color = {0, 0, 127}));
      connect(realExpression30.y, clark3p3.V_abc_D[3]) annotation(
        Line(points = {{88, 190}, {94, 190}, {94, 209}, {109, 209}}, color = {0, 0, 127}));
      connect(realExpression32.y, clark3p3.V_abc_D[1]) annotation(
        Line(points = {{88, 222}, {102, 222}, {102, 209}, {109, 209}}, color = {0, 0, 127}));
      connect(realExpression31.y, clark3p3.V_abc_D[2]) annotation(
        Line(points = {{88, 206}, {94, 206}, {94, 209}, {109, 209}}, color = {0, 0, 127}));
      connect(constant3.y, pcq.u2) annotation(
        Line(points = {{184, 186}, {184, 185}, {204, 185}, {204, 182}}, color = {0, 0, 127}));
      connect(pLL_fz21.beta, clark3p3.V_ab_D[2]) annotation(
        Line(points = {{161.8, 208.4}, {148.8, 208.4}, {148.8, 209.4}, {131.8, 209.4}}, color = {0, 0, 127}));
      connect(clark3p3.V_ab_D[1], pLL_fz21.alpha) annotation(
        Line(points = {{132, 209}, {144.5, 209}, {144.5, 217}, {162, 217}}, color = {0, 0, 127}));
      connect(pLL_fz21.q, pcq.u1) annotation(
        Line(points = {{185, 206.6}, {193, 206.6}, {193, 189.6}, {204, 189.6}}, color = {0, 0, 127}));
      connect(constant3.y, pcd.u2) annotation(
        Line(points = {{184, 186}, {197, 186}, {197, 206}, {205, 206}}, color = {0, 0, 127}));
      connect(pLL_fz21.d, pcd.u1) annotation(
        Line(points = {{185, 213}, {205, 213}, {205, 214}}, color = {0, 0, 127}));
      connect(ia.n, d_D_Transformer.pin) annotation(
        Line(points = {{242, 33}, {259, 33}, {259, -8}, {270, -8}}, color = {0, 0, 255}));
      connect(ib.n, d_D_Transformer.pin1) annotation(
        Line(points = {{244, -27}, {244, -26}, {270, -26}, {270, -21}}, color = {0, 0, 255}));
      connect(ic.n, d_D_Transformer.pin2) annotation(
        Line(points = {{246, -87}, {259, -87}, {259, -34}, {271, -34}}, color = {0, 0, 255}));
      connect(d_D_Transformer.pin3, ipca.p) annotation(
        Line(points = {{328, -7}, {333, -7}, {333, 38}, {345, 38}}, color = {0, 0, 255}));
      connect(pca.p, ipca.p) annotation(
        Line(points = {{368, 14}, {345, 14}, {345, 38}}, color = {0, 0, 255}));
      connect(ipcb.p, d_D_Transformer.pin4) annotation(
        Line(points = {{347, -22}, {337.5, -22}, {337.5, -21}, {328, -21}}, color = {0, 0, 255}));
      connect(pcb.p, ipcb.p) annotation(
        Line(points = {{371, -45}, {347, -45}, {347, -22}}, color = {0, 0, 255}));
      connect(ipcc.p, d_D_Transformer.pin5) annotation(
        Line(points = {{349, -82}, {334, -82}, {334, -34}, {327, -34}}, color = {0, 0, 255}));
      connect(pcc.p, ipcc.p) annotation(
        Line(points = {{366, -107}, {349, -107}, {349, -82}}, color = {0, 0, 255}));
      connect(realExpression47.y, division3.u1) annotation(
        Line(points = {{12, 287}, {50, 287}, {50, 289}}, color = {0, 0, 127}));
      connect(realExpression49.y, division2.u1) annotation(
        Line(points = {{12, 303}, {27, 303}, {27, 317}, {50, 317}}, color = {0, 0, 127}));
      connect(constant17.y, division4.u2) annotation(
        Line(points = {{32, 240}, {40, 240}, {40, 255}, {50, 255}}, color = {0, 0, 127}));
      connect(constant17.y, division3.u2) annotation(
        Line(points = {{32, 240}, {40, 240}, {40, 281}, {50, 281}}, color = {0, 0, 127}));
      connect(realExpression48.y, division4.u1) annotation(
        Line(points = {{12, 271}, {43, 271}, {43, 263}, {50, 263}}, color = {0, 0, 127}));
      connect(constant17.y, division2.u2) annotation(
        Line(points = {{32, 240}, {40, 240}, {40, 309}, {50, 309}}, color = {0, 0, 127}));
      connect(ramp.y, Controle_potencia.P_ref) annotation(
        Line(points = {{49, 171}, {81, 171}, {81, 113}, {89, 113}}, color = {0, 0, 127}));
      connect(ramp.y, Pref_pu.u1) annotation(
        Line(points = {{49, 171}, {115, 171}, {115, 163}, {122, 163}}, color = {0, 0, 127}));
      connect(ramp1.y, Controle_potencia.Q_ref) annotation(
        Line(points = {{86, 68}, {84, 68}, {84, 94}, {89, 94}}, color = {0, 0, 127}));
      connect(constant18.y, Qref_pu.u2) annotation(
        Line(points = {{125, 61}, {130, 61}, {130, 72}, {139, 72}}, color = {0, 0, 127}));
      connect(ramp1.y, Qref_pu.u1) annotation(
        Line(points = {{86, 68}, {96, 68}, {96, 80}, {139, 80}}, color = {0, 0, 127}));
      connect(constant12.y, pca_pu.u2) annotation(
        Line(points = {{605, 99}, {610, 99}, {610, 108}, {619, 108}}, color = {0, 0, 127}));
      connect(realExpression20.y, pca_pu.u1) annotation(
        Line(points = {{608, 119}, {613, 119}, {613, 116}, {619, 116}}, color = {0, 0, 127}));
      annotation(
        Diagram(coordinateSystem(extent = {{-150, 330}, {790, -210}}, grid = {1, 1})),
        Icon(coordinateSystem(extent = {{-120, -120}, {480, 160}}, grid = {1, 1})));
    end COBEP_fault;
    
    model SAEB1_S2_2
      //Base values
      parameter Real VAbase = 750e3;
      // RMS three phase
      parameter Real LVbase = 359.258;
      //phase peak value
      parameter Real MVbase = 11267.652;
      //phase peak value
      parameter Real ireflv = (VAbase/LVbase)*(2/3);
      //phase peak value
      parameter Real irefmv = (VAbase/MVbase)*(2/3);
      //phase peak value
      //BESS parameters
      parameter Real Vbat(start = 1024);
      //filter parameters
      parameter Real Rf = 1.5649e-4;
      parameter Real Lf = 1.6604e-5;
      parameter Real Cf = 0.0013;
      parameter Real Rd = 0.03;
      // grid parameters
      parameter Real Vth(start = 11.520e3);
      parameter Real Rr(start = 1.38003);
      parameter Real Lr(start = 0.003436);
      parameter Real f = 60.08;
      // Inputs
      parameter Real rampP(start = -750e3);
      parameter Real rampQ(start = 0);
      parameter Real rampP_duration(start = 0.432539);
      parameter Real rampQ_duration(start = 0);
      parameter Real rampP_start(start = 0.6);
      parameter Real rampQ_start(start = 0);
     //fault parameters
      parameter Real Rcc = 0.1;
      //Internal loads parameters
      parameter Real Rload = 9.207897;
      parameter Real Lload = 0.002754;
      
      VSC_FZ.Testes.VSC vsc annotation(
        Placement(visible = true, transformation(origin = {-138, 27}, extent = {{-13, -13}, {13, 13}}, rotation = 0)));
      Modelica.Electrical.Analog.Sources.ConstantVoltage constantVoltage(V = Vbat/2) annotation(
        Placement(visible = true, transformation(origin = {-223, 18}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      Modelica.Electrical.Analog.Basic.Resistor resistor(R = Rf, alpha = 0) annotation(
        Placement(visible = true, transformation(origin = {-50, 34}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Resistor resistor1(R = Rf, alpha = 0) annotation(
        Placement(visible = true, transformation(origin = {-47, -26}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Resistor resistor2(R = Rf, alpha = 0) annotation(
        Placement(visible = true, transformation(origin = {-45, -86}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Inductor inductor(L = Lf) annotation(
        Placement(visible = true, transformation(origin = {-22, 34}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Inductor inductor1(L = Lf) annotation(
        Placement(visible = true, transformation(origin = {-19, -26}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.VoltageSensor vsa annotation(
        Placement(visible = true, transformation(origin = {76, 52}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.VoltageSensor vsb annotation(
        Placement(visible = true, transformation(origin = {72, -5}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.VoltageSensor vscv annotation(
        Placement(visible = true, transformation(origin = {77, -61}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.CurrentSensor ia annotation(
        Placement(visible = true, transformation(origin = {99, 33}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.CurrentSensor ib annotation(
        Placement(visible = true, transformation(origin = {101, -27}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.CurrentSensor ic annotation(
        Placement(visible = true, transformation(origin = {103, -87}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
      VSC_FZ.Testes.Clark3p clark3p1 annotation(
        Placement(visible = true, transformation(origin = {229, 165}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression(y = ia.i) annotation(
        Placement(visible = true, transformation(origin = {195, 163}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression1(y = ib.i) annotation(
        Placement(visible = true, transformation(origin = {195, 147}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression2(y = ic.i) annotation(
        Placement(visible = true, transformation(origin = {195, 131}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression4(y = vsb.v) annotation(
        Placement(visible = true, transformation(origin = {-155, 159}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression5(y = vscv.v) annotation(
        Placement(visible = true, transformation(origin = {-155, 143}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression6(y = vsa.v) annotation(
        Placement(visible = true, transformation(origin = {-155, 175}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Sources.SineVoltage sineVoltage(V = Vth, freqHz = f, offset = 0, startTime = 0) annotation(
        Placement(visible = true, transformation(origin = {453, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Inductor inductor2(L = Lf) annotation(
        Placement(visible = true, transformation(origin = {-17, -86}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Sources.SineVoltage sineVoltage1(V = Vth, freqHz = f, offset = 0, phase = -2.094395102393195, startTime = 0) annotation(
        Placement(visible = true, transformation(origin = {455, -21}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Sources.SineVoltage sineVoltage2(V = Vth, freqHz = f, offset = 0, phase = 2.094395102393195, startTime = 0) annotation(
        Placement(visible = true, transformation(origin = {455, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Ground ground1 annotation(
        Placement(visible = true, transformation(origin = {490, -22}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
      Modelica.Blocks.Sources.RealExpression realExpression7(y = dsogi_pll.d) annotation(
        Placement(visible = true, transformation(origin = {368, 192.5}, extent = {{-17, -9.5}, {17, 9.5}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression8(y = dsogi_pll.q) annotation(
        Placement(visible = true, transformation(origin = {335.5, 117.5}, extent = {{-15.5, -9.5}, {15.5, 9.5}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression9(y = Controle_potencia.idref) annotation(
        Placement(visible = true, transformation(origin = {357, 177.5}, extent = {{-28, -9.5}, {28, 9.5}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression10(y = Correntes.d) annotation(
        Placement(visible = true, transformation(origin = {365, 164}, extent = {{-19, -10}, {19, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression11(y = Correntes.q) annotation(
        Placement(visible = true, transformation(origin = {365.5, 145}, extent = {{-19.5, -10}, {19.5, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression12(y = Controle_potencia.iqref) annotation(
        Placement(visible = true, transformation(origin = {359, 133.5}, extent = {{-26, -8.5}, {26, 8.5}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression13(y = dsogi_pll.ang) annotation(
        Placement(visible = true, transformation(origin = {479, 120.5}, extent = {{-23, -11.5}, {23, 11.5}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant2(k = 0) annotation(
        Placement(visible = true, transformation(origin = {523, 121}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression14(y = clarkInv.a) annotation(
        Placement(visible = true, transformation(origin = {-210, -56.5}, extent = {{-17, -8.5}, {17, 8.5}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression15(y = clarkInv.b) annotation(
        Placement(visible = true, transformation(origin = {-212, -83.5}, extent = {{-17, -9.5}, {17, 9.5}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression16(y = clarkInv.c) annotation(
        Placement(visible = true, transformation(origin = {-213, -111.5}, extent = {{-18, -9.5}, {18, 9.5}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.VoltageSensor vdc_2 annotation(
        Placement(visible = true, transformation(origin = {-246, 21}, extent = {{-7, -7}, {7, 7}}, rotation = -90)));
      Modelica.Blocks.Math.Division division annotation(
        Placement(visible = true, transformation(origin = {468, 181}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Math.Division division1 annotation(
        Placement(visible = true, transformation(origin = {466, 144}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Nonlinear.Limiter limiter(limitsAtInit = true, uMax = 1, uMin = -1) annotation(
        Placement(visible = true, transformation(origin = {-165, -58}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
      Modelica.Blocks.Nonlinear.Limiter limiter1(limitsAtInit = true, uMax = 1, uMin = -1) annotation(
        Placement(visible = true, transformation(origin = {-165, -84}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
      Modelica.Blocks.Nonlinear.Limiter limiter2(limitsAtInit = true, uMax = 1, uMin = -1) annotation(
        Placement(visible = true, transformation(origin = {-165, -112}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
      VSC_FZ.Testes.InversePark_fz inversePark_fz annotation(
        Placement(visible = true, transformation(origin = {508, 164}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      VSC_FZ.Testes.Park_Fz Correntes annotation(
        Placement(visible = true, transformation(origin = {267, 161}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Continuous.Filter filter(A_ripple = 0.1, analogFilter = Modelica.Blocks.Types.AnalogFilter.Bessel, f_cut = 60, filterType = Modelica.Blocks.Types.FilterType.LowPass) annotation(
        Placement(visible = true, transformation(origin = {595, -13}, extent = {{84, 4}, {104, 24}}, rotation = 0)));
      VSC_FZ.Testes.PowerCalc powerCalc annotation(
        Placement(visible = true, transformation(origin = {654, -6}, extent = {{-14, -14}, {14, 14}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression17(y = vsa.v) annotation(
        Placement(visible = true, transformation(origin = {620, 14}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression24(y = vsb.v) annotation(
        Placement(visible = true, transformation(origin = {620, 1}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression25(y = vscv.v) annotation(
        Placement(visible = true, transformation(origin = {620, -12}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression26(y = ia.i) annotation(
        Placement(visible = true, transformation(origin = {620, -27}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression27(y = ib.i) annotation(
        Placement(visible = true, transformation(origin = {620, -43}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression28(y = ic.i) annotation(
        Placement(visible = true, transformation(origin = {620, -59}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Continuous.Filter filter1(A_ripple = 0.1, analogFilter = Modelica.Blocks.Types.AnalogFilter.Bessel, f_cut = 60, filterType = Modelica.Blocks.Types.FilterType.LowPass) annotation(
        Placement(visible = true, transformation(origin = {595, -45}, extent = {{84, 4}, {104, 24}}, rotation = 0)));
      VSC_FZ.Testes.PowerControl1 Controle_potencia(Pctrl_KI = 0.168, Pctrl_KP = 5.36e-4) annotation(
        Placement(visible = true, transformation(origin = {96.5, 137.5}, extent = {{-17.5, -17.5}, {17.5, 17.5}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression23(y = filter.y) annotation(
        Placement(visible = true, transformation(origin = {21.5, 144}, extent = {{-20.5, -10}, {20.5, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression29(y = filter1.y) annotation(
        Placement(visible = true, transformation(origin = {21.5, 127}, extent = {{-20.5, -10}, {20.5, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression18(y = vdc_2.v) annotation(
        Placement(visible = true, transformation(origin = {428, 97.5}, extent = {{-23, -11.5}, {23, 11.5}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Resistor RG1(R = Rr, alpha = 0) annotation(
        Placement(visible = true, transformation(origin = {390, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Resistor RG2(R = Rr, alpha = 0) annotation(
        Placement(visible = true, transformation(origin = {389, -22}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Resistor RG3(R = Rr, alpha = 0) annotation(
        Placement(visible = true, transformation(origin = {389, -82}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Ground ground2 annotation(
        Placement(visible = true, transformation(origin = {92.5, -5.5}, extent = {{-6.5, -6.5}, {6.5, 6.5}}, rotation = 90)));
      Modelica.Electrical.Analog.Basic.Ground ground3 annotation(
        Placement(visible = true, transformation(origin = {94.5, 51.5}, extent = {{-6.5, -6.5}, {6.5, 6.5}}, rotation = 90)));
      Modelica.Electrical.Analog.Basic.Ground ground4 annotation(
        Placement(visible = true, transformation(origin = {99.5, -61.5}, extent = {{-6.5, -6.5}, {6.5, 6.5}}, rotation = 90)));
      Modelica.Electrical.Analog.Basic.Inductor LG1(L = Lr) annotation(
        Placement(visible = true, transformation(origin = {421, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Inductor LG2(L = Lr) annotation(
        Placement(visible = true, transformation(origin = {422, -21}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Inductor LG3(L = Lr) annotation(
        Placement(visible = true, transformation(origin = {423, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Ground ground5 annotation(
        Placement(visible = true, transformation(origin = {393.5, 13.5}, extent = {{-6.5, -6.5}, {6.5, 6.5}}, rotation = 90)));
      Modelica.Electrical.Analog.Sensors.VoltageSensor pca annotation(
        Placement(visible = true, transformation(origin = {375, 14}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.VoltageSensor pcc annotation(
        Placement(visible = true, transformation(origin = {373, -107}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Ground ground6 annotation(
        Placement(visible = true, transformation(origin = {391.5, -107.5}, extent = {{-6.5, -6.5}, {6.5, 6.5}}, rotation = 90)));
      Modelica.Electrical.Analog.Basic.Ground ground7 annotation(
        Placement(visible = true, transformation(origin = {396.5, -45.5}, extent = {{-6.5, -6.5}, {6.5, 6.5}}, rotation = 90)));
      Modelica.Electrical.Analog.Sensors.VoltageSensor pcb annotation(
        Placement(visible = true, transformation(origin = {378, -45}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Ground ground8 annotation(
        Placement(visible = true, transformation(origin = {-44.5, 9.5}, extent = {{-6.5, -6.5}, {6.5, 6.5}}, rotation = 90)));
      Modelica.Electrical.Analog.Sensors.VoltageSensor vta annotation(
        Placement(visible = true, transformation(origin = {-63, 10}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.VoltageSensor vtb annotation(
        Placement(visible = true, transformation(origin = {-59, -51}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Ground ground9 annotation(
        Placement(visible = true, transformation(origin = {-40.5, -51.5}, extent = {{-6.5, -6.5}, {6.5, 6.5}}, rotation = 90)));
      Modelica.Electrical.Analog.Basic.Ground ground10 annotation(
        Placement(visible = true, transformation(origin = {-35.5, -115.5}, extent = {{-6.5, -6.5}, {6.5, 6.5}}, rotation = 90)));
      Modelica.Electrical.Analog.Sensors.VoltageSensor vtc annotation(
        Placement(visible = true, transformation(origin = {-54, -115}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Math.Division vsq annotation(
        Placement(visible = true, transformation(origin = {-42, 117}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Math.Division vsd annotation(
        Placement(visible = true, transformation(origin = {-41, 141}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant4(k = LVbase) annotation(
        Placement(visible = true, transformation(origin = {-81, 116}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Division Pref_pu annotation(
        Placement(visible = true, transformation(origin = {118, 192}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant5(k = VAbase) annotation(
        Placement(visible = true, transformation(origin = {85, 177}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant6(k = VAbase) annotation(
        Placement(visible = true, transformation(origin = {695, 35}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Division P_pu annotation(
        Placement(visible = true, transformation(origin = {728, 51}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant7(k = VAbase) annotation(
        Placement(visible = true, transformation(origin = {695, -84}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Division Q_pu annotation(
        Placement(visible = true, transformation(origin = {728, -68}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant8(k = ireflv) annotation(
        Placement(visible = true, transformation(origin = {286, 194}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Division id_pu annotation(
        Placement(visible = true, transformation(origin = {320, 208}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant9(k = ireflv) annotation(
        Placement(visible = true, transformation(origin = {262, 114}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Division iq_pu annotation(
        Placement(visible = true, transformation(origin = {295, 130}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant10(k = ireflv) annotation(
        Placement(visible = true, transformation(origin = {879, -221}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Division idref_pu annotation(
        Placement(visible = true, transformation(origin = {912, -205}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression19(y = Controle_potencia.idref) annotation(
        Placement(visible = true, transformation(origin = {849, -193.5}, extent = {{-28, -9.5}, {28, 9.5}}, rotation = 0)));
      Modelica.Blocks.Math.Division iqref_pu annotation(
        Placement(visible = true, transformation(origin = {1026, -205}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant11(k = ireflv) annotation(
        Placement(visible = true, transformation(origin = {993, -221}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression21(y = Controle_potencia.iqref) annotation(
        Placement(visible = true, transformation(origin = {967, -196.5}, extent = {{-26, -8.5}, {26, 8.5}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.CurrentSensor ipca annotation(
        Placement(visible = true, transformation(origin = {351, 38}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.CurrentSensor ipcb annotation(
        Placement(visible = true, transformation(origin = {353, -22}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.CurrentSensor ipcc annotation(
        Placement(visible = true, transformation(origin = {355, -82}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant1(k = VAbase) annotation(
        Placement(visible = true, transformation(origin = {844, 39}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression22(y = ipcb.i) annotation(
        Placement(visible = true, transformation(origin = {769, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant13(k = VAbase) annotation(
        Placement(visible = true, transformation(origin = {844, -81}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Continuous.Filter filter2(analogFilter = Modelica.Blocks.Types.AnalogFilter.Bessel, f_cut = 60, filterType = Modelica.Blocks.Types.FilterType.LowPass) annotation(
        Placement(visible = true, transformation(origin = {744, -42}, extent = {{84, 4}, {104, 24}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression33(y = ipca.i) annotation(
        Placement(visible = true, transformation(origin = {769, -24}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Division P_pupc annotation(
        Placement(visible = true, transformation(origin = {877, 54}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Continuous.Filter filter3(analogFilter = Modelica.Blocks.Types.AnalogFilter.Bessel, f_cut = 60, filterType = Modelica.Blocks.Types.FilterType.LowPass) annotation(
        Placement(visible = true, transformation(origin = {744, -11}, extent = {{84, 4}, {104, 24}}, rotation = 0)));
      Modelica.Blocks.Math.Division Q_pupc annotation(
        Placement(visible = true, transformation(origin = {877, -65}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression34(y = pcb.v) annotation(
        Placement(visible = true, transformation(origin = {769, 4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression35(y = ipcc.i) annotation(
        Placement(visible = true, transformation(origin = {769, -56}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression36(y = pcc.v) annotation(
        Placement(visible = true, transformation(origin = {769, -9}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      VSC_FZ.Testes.PowerCalc powerCalc_pcc annotation(
        Placement(visible = true, transformation(origin = {803, -4}, extent = {{-14, -14}, {14, 14}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression37(y = pca.v) annotation(
        Placement(visible = true, transformation(origin = {769, 17}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Ideal.IdealClosingSwitch idealClosingSwitch annotation(
        Placement(visible = true, transformation(origin = {428, -48}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
      Modelica.Electrical.Analog.Basic.Resistor resistor6(R = Rcc, alpha = 0) annotation(
        Placement(visible = true, transformation(origin = {443, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Ground ground11 annotation(
        Placement(visible = true, transformation(origin = {463.5, 0.5}, extent = {{-6.5, -6.5}, {6.5, 6.5}}, rotation = 90)));
      Modelica.Electrical.Analog.Ideal.IdealClosingSwitch idealClosingSwitch1 annotation(
        Placement(visible = true, transformation(origin = {428, -106}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
      Modelica.Electrical.Analog.Ideal.IdealClosingSwitch switch2 annotation(
        Placement(visible = true, transformation(origin = {425, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
      Modelica.Blocks.MathBoolean.And and1(nu = 2) annotation(
        Placement(visible = true, transformation(origin = {382, -184}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
      Modelica.Blocks.Sources.BooleanStep booleanStep4(startTime = 0, startValue = true) annotation(
        Placement(visible = true, transformation(origin = {355, -169}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Ground ground13 annotation(
        Placement(visible = true, transformation(origin = {465.5, -115.5}, extent = {{-6.5, -6.5}, {6.5, 6.5}}, rotation = 90)));
      Modelica.Electrical.Analog.Basic.Resistor resistor7(R = Rcc, alpha = 0) annotation(
        Placement(visible = true, transformation(origin = {445, -58}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Ground ground12 annotation(
        Placement(visible = true, transformation(origin = {465.5, -57.5}, extent = {{-6.5, -6.5}, {6.5, 6.5}}, rotation = 90)));
      Modelica.Blocks.Sources.BooleanStep booleanStep3(startTime = 0, startValue = false) annotation(
        Placement(visible = true, transformation(origin = {354, -202}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Resistor resistor8(R = Rcc, alpha = 0) annotation(
        Placement(visible = true, transformation(origin = {445, -118}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.BooleanStep booleanStep(startTime = 0, startValue = true) annotation(
        Placement(visible = true, transformation(origin = {282, -169}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.MathBoolean.And and2(nu = 2) annotation(
        Placement(visible = true, transformation(origin = {308, -184}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
      Modelica.Blocks.Sources.BooleanStep booleanStep1(startTime = 0, startValue = false) annotation(
        Placement(visible = true, transformation(origin = {280, -202}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.BooleanStep booleanStep2(startTime = 0, startValue = true) annotation(
        Placement(visible = true, transformation(origin = {401, -214}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.MathBoolean.And and3(nu = 2) annotation(
        Placement(visible = true, transformation(origin = {430, -201}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
      Modelica.Blocks.Sources.BooleanStep booleanStep5(startTime = 0, startValue = false) annotation(
        Placement(visible = true, transformation(origin = {399, -247}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Division ib_pu annotation(
        Placement(visible = true, transformation(origin = {788, -154}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression46(y = ia.i) annotation(
        Placement(visible = true, transformation(origin = {731, -136}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression44(y = ib.i) annotation(
        Placement(visible = true, transformation(origin = {731, -152}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression45(y = ic.i) annotation(
        Placement(visible = true, transformation(origin = {731, -168}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant16(k = ireflv) annotation(
        Placement(visible = true, transformation(origin = {751, -199}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Division ia_pu annotation(
        Placement(visible = true, transformation(origin = {788, -127}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Math.Division ic_pu annotation(
        Placement(visible = true, transformation(origin = {788, -180}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      VSC_FZ.Testes.ClarkInv clarkInv annotation(
        Placement(visible = true, transformation(origin = {574, 169}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Capacitor capacitor(C = Cf) annotation(
        Placement(visible = true, transformation(origin = {-4, 20}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      Modelica.Electrical.Analog.Basic.Resistor resistor9(R = Rd, alpha = 0) annotation(
        Placement(visible = true, transformation(origin = {17, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Resistor resistor10(R = Rf, alpha = 0) annotation(
        Placement(visible = true, transformation(origin = {21, 34}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Inductor inductor6(L = Lf) annotation(
        Placement(visible = true, transformation(origin = {49.5, 33.5}, extent = {{-10.5, -10.5}, {10.5, 10.5}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Resistor resistor11(R = Rf, alpha = 0) annotation(
        Placement(visible = true, transformation(origin = {22, -25}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Capacitor capacitor1(C = Cf) annotation(
        Placement(visible = true, transformation(origin = {-2, -40}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      Modelica.Electrical.Analog.Basic.Inductor inductor7(L = Lf) annotation(
        Placement(visible = true, transformation(origin = {50, -25}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Resistor resistor12(R = Rd, alpha = 0) annotation(
        Placement(visible = true, transformation(origin = {18, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Resistor resistor13(R = Rf, alpha = 0) annotation(
        Placement(visible = true, transformation(origin = {25, -86}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Capacitor capacitor2(C = Cf) annotation(
        Placement(visible = true, transformation(origin = {0, -100}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      Modelica.Electrical.Analog.Basic.Inductor inductor8(L = Lf) annotation(
        Placement(visible = true, transformation(origin = {53, -86}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Resistor resistor14(R = Rd, alpha = 0) annotation(
        Placement(visible = true, transformation(origin = {21, -110}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      VSC_FZ.Testes.DSOGI_PLL dsogi_pll annotation(
        Placement(visible = true, transformation(origin = {-107, 160}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      VSC_FZ.Testes.InnerControl_PIR innerControl_PIR(Ictrl_KI = 0.4916, Ictrl_KP = 5.216e-2, KIh = 600, freq = 2*3.1416*60, h = 2, lw = 2*3.1416*60*2*Lf, r = 2*Rf) annotation(
        Placement(visible = true, transformation(origin = {425.5, 169.5}, extent = {{-15.5, -15.5}, {15.5, 15.5}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression3(y = dsogi_pll.ang) annotation(
        Placement(visible = true, transformation(origin = {235.5, 137}, extent = {{-20.5, -9}, {20.5, 9}}, rotation = 0)));
      VSC_FZ.Testes.D_D_Transformer d_D_Transformer(L_trans_prim = 32.86e-3, L_trans_secon = 3.34e-5, Lm1 = 0, R_trans_prim = 1.1446, R_trans_secon = 11.638e-4, Vac_primary = 13800, Vac_secondary = 440, considerMagnetization = false, trans_ratio = 31.363636) annotation(
        Placement(visible = true, transformation(origin = {299, -30}, extent = {{27, -27}, {-27, 27}}, rotation = 0)));
      Modelica.Blocks.Math.Division ipca_pu annotation(
        Placement(visible = true, transformation(origin = {691, -135}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Math.Division ipcb_pu annotation(
        Placement(visible = true, transformation(origin = {691, -163}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression47(y = ipcb.i) annotation(
        Placement(visible = true, transformation(origin = {634, -161}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Division ipcc_pu annotation(
        Placement(visible = true, transformation(origin = {691, -189}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression48(y = ipcc.i) annotation(
        Placement(visible = true, transformation(origin = {634, -177}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression49(y = ipca.i) annotation(
        Placement(visible = true, transformation(origin = {634, -145}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant17(k = irefmv) annotation(
        Placement(visible = true, transformation(origin = {654, -208}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Ramp Pramp(duration = rampP_duration, height = rampP, startTime = rampP_start) annotation(
        Placement(visible = true, transformation(origin = {22, 206}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Ramp Qramp(duration = rampQ_duration, height = rampQ, startTime = rampQ_start) annotation(
        Placement(visible = true, transformation(origin = {32, 102}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant18(k = VAbase) annotation(
        Placement(visible = true, transformation(origin = {93, 90}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Division Qref_pu annotation(
        Placement(visible = true, transformation(origin = {135, 109}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Math.Division pca_pu annotation(
        Placement(visible = true, transformation(origin = {884, -138}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant12(k = MVbase) annotation(
        Placement(visible = true, transformation(origin = {851, -151}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression20(y = pca.v) annotation(
        Placement(visible = true, transformation(origin = {854, -131}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Sources.ConstantVoltage constantVoltage1(V = Vbat/2) annotation(
        Placement(visible = true, transformation(origin = {-223, -16}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      Modelica.Electrical.Analog.Sensors.CurrentSensor currentSensor annotation(
        Placement(visible = true, transformation(origin = {-234, 1}, extent = {{-4, -4}, {4, 4}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Ground ground annotation(
        Placement(visible = true, transformation(origin = {-202, 70}, extent = {{-54, -86}, {-38, -70}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.CurrentSensor filtro_a annotation(
        Placement(visible = true, transformation(origin = {39, 10}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.CurrentSensor filtro_b annotation(
        Placement(visible = true, transformation(origin = {39, -50}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.CurrentSensor filtro_c annotation(
        Placement(visible = true, transformation(origin = {42, -110}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.CurrentSensor filtro_terra annotation(
        Placement(visible = true, transformation(origin = {63, -111}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Ground ground14 annotation(
        Placement(visible = true, transformation(origin = {82.5, -111.5}, extent = {{-6.5, -6.5}, {6.5, 6.5}}, rotation = 90)));
      Modelica.Blocks.Math.Division vsa_pu annotation(
        Placement(visible = true, transformation(origin = {954, -138}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant constant19(k = LVbase) annotation(
        Placement(visible = true, transformation(origin = {921, -151}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression50(y = vsa.v) annotation(
        Placement(visible = true, transformation(origin = {924, -131}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.CurrentSensor Rloadb annotation(
        Placement(visible = true, transformation(origin = {200, -40}, extent = {{-6, -6}, {6, 6}}, rotation = -90)));
      Modelica.Electrical.Analog.Sensors.CurrentSensor Loadgb annotation(
        Placement(visible = true, transformation(origin = {200, -79}, extent = {{-6, -6}, {6, 6}}, rotation = -90)));
      Modelica.Electrical.Analog.Sensors.CurrentSensor Rloadc annotation(
        Placement(visible = true, transformation(origin = {145, -37}, extent = {{-6, -6}, {6, 6}}, rotation = -90)));
      Modelica.Electrical.Analog.Sensors.CurrentSensor Loadgc annotation(
        Placement(visible = true, transformation(origin = {145, -77}, extent = {{-6, -6}, {6, 6}}, rotation = -90)));
      Modelica.Electrical.Analog.Sensors.VoltageSensor vsab annotation(
        Placement(visible = true, transformation(origin = {115, 7}, extent = {{-7, -7}, {7, 7}}, rotation = -90)));
      Modelica.Electrical.Analog.Sensors.VoltageSensor vsbc annotation(
        Placement(visible = true, transformation(origin = {116, -59}, extent = {{-7, -7}, {7, 7}}, rotation = -90)));
      Modelica.Electrical.Analog.Sensors.VoltageSensor vsac annotation(
        Placement(visible = true, transformation(origin = {123, -16}, extent = {{-7, -7}, {7, 7}}, rotation = -90)));
    Modelica.Electrical.Analog.Basic.Resistor resistor19(R = Rload, alpha = 0) annotation(
        Placement(visible = true, transformation(origin = {223, -64}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Electrical.Analog.Basic.Resistor resistor20(R = Rload, alpha = 0) annotation(
        Placement(visible = true, transformation(origin = {166, -72}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Blocks.Math.Division pcb_pu annotation(
        Placement(visible = true, transformation(origin = {941, -79}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
    Modelica.Blocks.Math.Division vsb_pu annotation(
        Placement(visible = true, transformation(origin = {1011, -79}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
    Modelica.Blocks.Sources.Constant constant3(k = MVbase) annotation(
        Placement(visible = true, transformation(origin = {908, -92}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.RealExpression realExpression30(y = vsb.v) annotation(
        Placement(visible = true, transformation(origin = {981, -72}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.RealExpression realExpression31(y = pcb.v) annotation(
        Placement(visible = true, transformation(origin = {911, -72}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Constant constant14(k = LVbase) annotation(
        Placement(visible = true, transformation(origin = {978, -92}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.RealExpression realExpression32(y = pcc.v) annotation(
        Placement(visible = true, transformation(origin = {756, -260}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Constant constant15(k = LVbase) annotation(
        Placement(visible = true, transformation(origin = {823, -280}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Math.Division pcc_pu annotation(
        Placement(visible = true, transformation(origin = {786, -267}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
    Modelica.Blocks.Sources.Constant constant20(k = MVbase) annotation(
        Placement(visible = true, transformation(origin = {753, -280}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.RealExpression realExpression38(y = vscv.v) annotation(
        Placement(visible = true, transformation(origin = {826, -260}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Math.Division vscv_pu annotation(
        Placement(visible = true, transformation(origin = {856, -267}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
    Modelica.Electrical.Analog.Sensors.CurrentSensor currentSensor1 annotation(
        Placement(visible = true, transformation(origin = {196, 25}, extent = {{-6, -6}, {6, 6}}, rotation = -90)));
    Modelica.Electrical.Analog.Basic.Resistor resistor15(R = Rload, alpha = 0) annotation(
        Placement(visible = true, transformation(origin = {217, -8}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Electrical.Analog.Sensors.CurrentSensor currentSensor2 annotation(
        Placement(visible = true, transformation(origin = {196, -15}, extent = {{-6, -6}, {6, 6}}, rotation = -90)));
    Modelica.Electrical.Analog.Basic.Inductor inductor3(L = Lload, i(start = 0)) annotation(
        Placement(visible = true, transformation(origin = {217, 15}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Electrical.Analog.Basic.Inductor inductor4(L = Lload, i(start = 0)) annotation(
        Placement(visible = true, transformation(origin = {223, -40}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Electrical.Analog.Basic.Inductor inductor5(L = Lload, i(start = 0)) annotation(
        Placement(visible = true, transformation(origin = {166, -48}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    equation
      connect(constantVoltage.p, vsc.pin_p) annotation(
        Line(points = {{-223, 28}, {-223, 34}, {-152, 34}, {-152, 33.5}}, color = {0, 0, 255}));
      connect(resistor.n, inductor.p) annotation(
        Line(points = {{-40, 34}, {-32, 34}}, color = {0, 0, 255}));
      connect(resistor1.n, inductor1.p) annotation(
        Line(points = {{-37, -26}, {-29, -26}}, color = {0, 0, 255}));
      connect(realExpression.y, clark3p1.V_abc_D[1]) annotation(
        Line(points = {{206, 163}, {211, 163}, {211, 165}, {217, 165}}, color = {0, 0, 127}));
      connect(realExpression1.y, clark3p1.V_abc_D[2]) annotation(
        Line(points = {{206, 147}, {208, 147}, {208, 165}, {217, 165}}, color = {0, 0, 127}));
      connect(realExpression2.y, clark3p1.V_abc_D[3]) annotation(
        Line(points = {{206, 131}, {208, 131}, {208, 165}, {217, 165}}, color = {0, 0, 127}));
      connect(inductor2.p, resistor2.n) annotation(
        Line(points = {{-27, -86}, {-35, -86}}, color = {0, 0, 255}));
      connect(sineVoltage.n, ground1.p) annotation(
        Line(points = {{463, 38}, {480, 38}, {480, -22}}, color = {0, 0, 255}));
      connect(sineVoltage1.n, ground1.p) annotation(
        Line(points = {{465, -21}, {472.5, -21}, {472.5, -22}, {480, -22}}, color = {0, 0, 255}));
      connect(sineVoltage2.n, ground1.p) annotation(
        Line(points = {{465, -80}, {465, -82}, {480, -82}, {480, -22}}, color = {0, 0, 255}));
      connect(vdc_2.p, constantVoltage.p) annotation(
        Line(points = {{-246, 28}, {-224, 28}}, color = {0, 0, 255}));
      connect(vdc_2.n, constantVoltage.n) annotation(
        Line(points = {{-246, 14}, {-246, 8}, {-224, 8}}, color = {0, 0, 255}));
      connect(realExpression14.y, limiter.u) annotation(
        Line(points = {{-191.3, -56.5}, {-183.3, -56.5}, {-183.3, -58}, {-175.3, -58}}, color = {0, 0, 127}));
      connect(limiter.y, vsc.ma) annotation(
        Line(points = {{-156.2, -58}, {-145.2, -58}, {-145.2, 12}}, color = {0, 0, 127}));
      connect(realExpression15.y, limiter1.u) annotation(
        Line(points = {{-193.3, -83.5}, {-184.3, -83.5}, {-184.3, -84}, {-175.3, -84}}, color = {0, 0, 127}));
      connect(realExpression16.y, limiter2.u) annotation(
        Line(points = {{-193.2, -111.5}, {-184.2, -111.5}, {-184.2, -112}, {-175.2, -112}}, color = {0, 0, 127}));
      connect(limiter1.y, vsc.mb) annotation(
        Line(points = {{-156.2, -84}, {-136.2, -84}, {-136.2, 12}}, color = {0, 0, 127}));
      connect(limiter2.y, vsc.mc) annotation(
        Line(points = {{-156.2, -112}, {-130.2, -112}, {-130.2, 12}}, color = {0, 0, 127}));
      connect(division.y, inversePark_fz.d) annotation(
        Line(points = {{475.7, 181}, {487.4, 181}, {487.4, 168}, {495.4, 168}}, color = {0, 0, 127}));
      connect(division1.y, inversePark_fz.q) annotation(
        Line(points = {{473.7, 144}, {489.7, 144}, {489.7, 160}, {495.7, 160}}, color = {0, 0, 127}));
      connect(inversePark_fz.theta, realExpression13.y) annotation(
        Line(points = {{508, 152}, {508, 120.5}, {504, 120.5}}, color = {0, 0, 127}));
      connect(clark3p1.V_ab_D[1], Correntes.alpha) annotation(
        Line(points = {{240, 165}, {255, 165}}, color = {0, 0, 127}));
      connect(clark3p1.V_ab_D[2], Correntes.beta) annotation(
        Line(points = {{240, 165}, {248, 165}, {248, 157}, {255, 157}}, color = {0, 0, 127}));
      connect(realExpression17.y, powerCalc.V1) annotation(
        Line(points = {{631, 14}, {636, 14}, {636, 6}, {639, 6}}, color = {0, 0, 127}));
      connect(realExpression24.y, powerCalc.V2) annotation(
        Line(points = {{631, 1}, {635, 1}, {635, 2}, {639, 2}}, color = {0, 0, 127}));
      connect(realExpression25.y, powerCalc.V3) annotation(
        Line(points = {{631, -12}, {633, -12}, {633, -2}, {639, -2}}, color = {0, 0, 127}));
      connect(realExpression26.y, powerCalc.I1) annotation(
        Line(points = {{631, -27}, {635, -27}, {635, -9}, {639, -9}}, color = {0, 0, 127}));
      connect(realExpression27.y, powerCalc.I2) annotation(
        Line(points = {{631, -43}, {636, -43}, {636, -13}, {639, -13}}, color = {0, 0, 127}));
      connect(realExpression28.y, powerCalc.I3) annotation(
        Line(points = {{631, -59}, {637, -59}, {637, -18}, {639, -18}}, color = {0, 0, 127}));
      connect(powerCalc.P, filter.u) annotation(
        Line(points = {{669.4, 0.44}, {669.4, 1.44}, {677.4, 1.44}}, color = {0, 0, 127}));
      connect(powerCalc.Q, filter1.u) annotation(
        Line(points = {{669.4, -10.2}, {672.8, -10.2}, {672.8, -31.4}, {677.8, -31.4}}, color = {0, 0, 127}));
      connect(realExpression29.y, Controle_potencia.Q) annotation(
        Line(points = {{44.05, 127}, {60.05, 127}, {60.05, 133}, {77.05, 133}}, color = {0, 0, 127}));
      connect(realExpression23.y, Controle_potencia.P) annotation(
        Line(points = {{44.05, 144}, {55.05, 144}, {55.05, 139}, {77.05, 139}}, color = {0, 0, 127}));
      connect(realExpression18.y, division1.u2) annotation(
        Line(points = {{453.3, 97.5}, {453.3, 140}, {458.3, 140}}, color = {0, 0, 127}));
      connect(realExpression18.y, division.u2) annotation(
        Line(points = {{453.3, 97.5}, {453.3, 177}, {460.3, 177}}, color = {0, 0, 127}));
      connect(vsb.n, ground2.p) annotation(
        Line(points = {{79, -5}, {86, -5}}, color = {0, 0, 255}));
      connect(vsa.n, ground3.p) annotation(
        Line(points = {{83, 52}, {88, 52}}, color = {0, 0, 255}));
      connect(vscv.n, ground4.p) annotation(
        Line(points = {{85, -61}, {93, -61}}, color = {0, 0, 255}));
      connect(RG1.n, LG1.p) annotation(
        Line(points = {{400, 38}, {404, 38}, {404, 40}, {411, 40}}, color = {0, 0, 255}));
      connect(LG1.n, sineVoltage.p) annotation(
        Line(points = {{431, 40}, {437, 40}, {437, 38}, {443, 38}}, color = {0, 0, 255}));
      connect(pca.n, ground5.p) annotation(
        Line(points = {{382, 14}, {387, 14}}, color = {0, 0, 255}));
      connect(pcc.n, ground6.p) annotation(
        Line(points = {{380, -107}, {385, -107}}, color = {0, 0, 255}));
      connect(pcb.n, ground7.p) annotation(
        Line(points = {{385, -45}, {390, -45}}, color = {0, 0, 255}));
      connect(vta.n, ground8.p) annotation(
        Line(points = {{-56, 10}, {-51, 10}}, color = {0, 0, 255}));
      connect(vtb.n, ground9.p) annotation(
        Line(points = {{-52, -51}, {-47, -51}}, color = {0, 0, 255}));
      connect(vtc.n, ground10.p) annotation(
        Line(points = {{-47, -115}, {-42, -115}}, color = {0, 0, 255}));
      connect(vtc.p, resistor2.p) annotation(
        Line(points = {{-61, -115}, {-67, -115}, {-67, -86}, {-55, -86}}, color = {0, 0, 255}));
      connect(vtb.p, resistor1.p) annotation(
        Line(points = {{-66, -51}, {-72, -51}, {-72, -26}, {-57, -26}}, color = {0, 0, 255}));
      connect(vta.p, resistor.p) annotation(
        Line(points = {{-70, 10}, {-73, 10}, {-73, 34}, {-60, 34}}, color = {0, 0, 255}));
      connect(constant4.y, vsq.u2) annotation(
        Line(points = {{-70, 116}, {-50, 116}, {-50, 113}}, color = {0, 0, 127}));
      connect(constant4.y, vsd.u2) annotation(
        Line(points = {{-70, 116}, {-57, 116}, {-57, 137}, {-49, 137}}, color = {0, 0, 127}));
      connect(constant5.y, Pref_pu.u2) annotation(
        Line(points = {{96, 177}, {101, 177}, {101, 188}, {110, 188}}, color = {0, 0, 127}));
      connect(constant6.y, P_pu.u2) annotation(
        Line(points = {{706, 35}, {711, 35}, {711, 47}, {720, 47}}, color = {0, 0, 127}));
      connect(powerCalc.P, P_pu.u1) annotation(
        Line(points = {{669.4, 0.44}, {671.8, 0.44}, {671.8, 55.88}, {720.8, 55.88}}, color = {0, 0, 127}));
      connect(constant7.y, Q_pu.u2) annotation(
        Line(points = {{706, -84}, {711, -84}, {711, -72}, {720, -72}}, color = {0, 0, 127}));
      connect(powerCalc.Q, Q_pu.u1) annotation(
        Line(points = {{669.4, -10.2}, {673.8, -10.2}, {673.8, -64.4}, {720.8, -64.4}}, color = {0, 0, 127}));
      connect(constant8.y, id_pu.u2) annotation(
        Line(points = {{297, 194}, {297, 204}, {312, 204}}, color = {0, 0, 127}));
      connect(id_pu.u1, Correntes.d) annotation(
        Line(points = {{311.6, 212.2}, {265.6, 212.2}, {265.6, 171.2}, {286.6, 171.2}, {286.6, 165}, {278, 165}}, color = {0, 0, 127}));
      connect(constant9.y, iq_pu.u2) annotation(
        Line(points = {{273, 114}, {278, 114}, {278, 126}, {287, 126}}, color = {0, 0, 127}));
      connect(iq_pu.u1, Correntes.q) annotation(
        Line(points = {{286.6, 134.2}, {280.6, 134.2}, {280.6, 157}, {278, 157}}, color = {0, 0, 127}));
      connect(constant10.y, idref_pu.u2) annotation(
        Line(points = {{890, -221}, {895, -221}, {895, -209}, {904, -209}}, color = {0, 0, 127}));
      connect(realExpression19.y, idref_pu.u1) annotation(
        Line(points = {{879.8, -193.5}, {897.8, -193.5}, {897.8, -201.5}, {903.8, -201.5}}, color = {0, 0, 127}));
      connect(constant11.y, iqref_pu.u2) annotation(
        Line(points = {{1004, -221}, {1009, -221}, {1009, -209}, {1018, -209}}, color = {0, 0, 127}));
      connect(realExpression21.y, iqref_pu.u1) annotation(
        Line(points = {{995.6, -196.5}, {1006.6, -196.5}, {1006.6, -201}, {1017.6, -201}}, color = {0, 0, 127}));
      connect(RG1.p, ipca.n) annotation(
        Line(points = {{380, 38}, {357, 38}}, color = {0, 0, 255}));
      connect(RG2.p, ipcb.n) annotation(
        Line(points = {{379, -22}, {359, -22}}, color = {0, 0, 255}));
      connect(RG3.p, ipcc.n) annotation(
        Line(points = {{379, -82}, {361, -82}}, color = {0, 0, 255}));
      connect(constant1.y, P_pupc.u2) annotation(
        Line(points = {{855, 39}, {860, 39}, {860, 50}, {869, 50}}, color = {0, 0, 127}));
      connect(realExpression33.y, powerCalc_pcc.I1) annotation(
        Line(points = {{780, -24}, {784, -24}, {784, -7}, {788, -7}}, color = {0, 0, 127}));
      connect(realExpression35.y, powerCalc_pcc.I3) annotation(
        Line(points = {{780, -56}, {786, -56}, {786, -16}, {788, -16}}, color = {0, 0, 127}));
      connect(constant13.y, Q_pupc.u2) annotation(
        Line(points = {{855, -81}, {860, -81}, {860, -69}, {869, -69}}, color = {0, 0, 127}));
      connect(realExpression37.y, powerCalc_pcc.V1) annotation(
        Line(points = {{780, 17}, {785, 17}, {785, 8}, {788, 8}}, color = {0, 0, 127}));
      connect(powerCalc_pcc.Q, Q_pupc.u1) annotation(
        Line(points = {{818.4, -8.2}, {822.4, -8.2}, {822.4, -61.2}, {869.4, -61.2}}, color = {0, 0, 127}));
      connect(realExpression22.y, powerCalc_pcc.I2) annotation(
        Line(points = {{780, -40}, {785, -40}, {785, -11}, {788, -11}}, color = {0, 0, 127}));
      connect(realExpression34.y, powerCalc_pcc.V2) annotation(
        Line(points = {{780, 4}, {788, 4}}, color = {0, 0, 127}));
      connect(realExpression36.y, powerCalc_pcc.V3) annotation(
        Line(points = {{780, -9}, {782, -9}, {782, 0}, {788, 0}}, color = {0, 0, 127}));
      connect(powerCalc_pcc.Q, filter2.u) annotation(
        Line(points = {{818.4, -8.2}, {821.4, -8.2}, {821.4, -28.2}, {826.4, -28.2}}, color = {0, 0, 127}));
      connect(powerCalc_pcc.P, filter3.u) annotation(
        Line(points = {{818.4, 2.44}, {826.4, 2.44}, {826.4, 3.44}}, color = {0, 0, 127}));
      connect(powerCalc_pcc.P, P_pupc.u1) annotation(
        Line(points = {{818.4, 2.44}, {820.4, 2.44}, {820.4, 58.44}, {869.4, 58.44}}, color = {0, 0, 127}));
      connect(booleanStep3.y, and1.u[2]) annotation(
        Line(points = {{365, -202}, {373, -202}, {373, -184}, {376, -184}}, color = {255, 0, 255}));
      connect(booleanStep4.y, and1.u[1]) annotation(
        Line(points = {{366, -169}, {373, -169}, {373, -184}, {376, -184}}, color = {255, 0, 255}));
      connect(switch2.n, RG1.p) annotation(
        Line(points = {{425, 20}, {425, 25}, {380, 25}, {380, 38}}, color = {0, 0, 255}));
      connect(resistor6.p, switch2.p) annotation(
        Line(points = {{433, 0}, {425, 0}}, color = {0, 0, 255}));
      connect(ground11.p, resistor6.n) annotation(
        Line(points = {{457, 0.5}, {453, 0.5}, {453, -0.5}}, color = {0, 0, 255}));
      connect(idealClosingSwitch.n, RG2.p) annotation(
        Line(points = {{428, -38}, {379, -38}, {379, -22}}, color = {0, 0, 255}));
      connect(idealClosingSwitch.p, resistor7.p) annotation(
        Line(points = {{428, -58}, {435, -58}}, color = {0, 0, 255}));
      connect(resistor7.n, ground12.p) annotation(
        Line(points = {{455, -58}, {459, -58}, {459, -57}}, color = {0, 0, 255}));
      connect(idealClosingSwitch1.n, RG3.p) annotation(
        Line(points = {{428, -96}, {375, -96}, {375, -82}, {379, -82}}, color = {0, 0, 255}));
      connect(idealClosingSwitch1.p, resistor8.p) annotation(
        Line(points = {{428, -116}, {435, -116}, {435, -118}}, color = {0, 0, 255}));
      connect(resistor8.n, ground13.p) annotation(
        Line(points = {{455, -118}, {459, -118}, {459, -115}}, color = {0, 0, 255}));
      connect(booleanStep.y, and2.u[1]) annotation(
        Line(points = {{293, -169}, {299, -169}, {299, -184}, {302, -184}}, color = {255, 0, 255}));
      connect(booleanStep1.y, and2.u[2]) annotation(
        Line(points = {{291, -202}, {299, -202}, {299, -184}, {302, -184}}, color = {255, 0, 255}));
      connect(booleanStep2.y, and3.u[1]) annotation(
        Line(points = {{412, -214}, {421, -214}, {421, -201}, {424, -201}}, color = {255, 0, 255}));
      connect(booleanStep5.y, and3.u[2]) annotation(
        Line(points = {{410, -247}, {421, -247}, {421, -201}, {424, -201}}, color = {255, 0, 255}));
      connect(constant16.y, ib_pu.u2) annotation(
        Line(points = {{762, -199}, {770, -199}, {770, -158}, {780, -158}}, color = {0, 0, 127}));
      connect(constant16.y, ia_pu.u2) annotation(
        Line(points = {{762, -199}, {770, -199}, {770, -131}, {780, -131}}, color = {0, 0, 127}));
      connect(realExpression45.y, ic_pu.u1) annotation(
        Line(points = {{742, -168}, {773, -168}, {773, -176}, {780, -176}}, color = {0, 0, 127}));
      connect(realExpression46.y, ia_pu.u1) annotation(
        Line(points = {{742, -136}, {757, -136}, {757, -123}, {780, -123}}, color = {0, 0, 127}));
      connect(constant16.y, ic_pu.u2) annotation(
        Line(points = {{762, -199}, {770, -199}, {770, -184}, {780, -184}}, color = {0, 0, 127}));
      connect(realExpression44.y, ib_pu.u1) annotation(
        Line(points = {{742, -152}, {780, -152}, {780, -150}}, color = {0, 0, 127}));
      connect(inversePark_fz.alpha, clarkInv.A) annotation(
        Line(points = {{519, 168}, {540, 168}, {540, 174}, {561, 174}}, color = {0, 0, 127}));
      connect(inversePark_fz.beta, clarkInv.B) annotation(
        Line(points = {{519, 160}, {548, 160}, {548, 168}, {561, 168}}, color = {0, 0, 127}));
      connect(constant2.y, clarkInv.C) annotation(
        Line(points = {{534, 121}, {554, 121}, {554, 163}, {561, 163}}, color = {0, 0, 127}));
      connect(capacitor.n, resistor9.p) annotation(
        Line(points = {{-4, 10}, {7, 10}}, color = {0, 0, 255}));
      connect(capacitor.p, inductor.n) annotation(
        Line(points = {{-4, 30}, {-4, 34}, {-12, 34}}, color = {0, 0, 255}));
      connect(capacitor.p, resistor10.p) annotation(
        Line(points = {{-4, 30}, {-4, 34}, {11, 34}}, color = {0, 0, 255}));
      connect(resistor10.n, inductor6.p) annotation(
        Line(points = {{31, 34}, {35, 34}, {35, 33.5}, {39, 33.5}}, color = {0, 0, 255}));
      connect(inductor1.n, resistor11.p) annotation(
        Line(points = {{-9, -26}, {12, -26}, {12, -25}}, color = {0, 0, 255}));
      connect(capacitor1.p, inductor1.n) annotation(
        Line(points = {{-2, -30}, {-2, -26}, {-9, -26}}, color = {0, 0, 255}));
      connect(capacitor1.n, resistor12.p) annotation(
        Line(points = {{-2, -50}, {8, -50}}, color = {0, 0, 255}));
      connect(resistor11.n, inductor7.p) annotation(
        Line(points = {{32, -25}, {40, -25}}, color = {0, 0, 255}));
      connect(inductor2.n, resistor13.p) annotation(
        Line(points = {{-7, -86}, {15, -86}}, color = {0, 0, 255}));
      connect(capacitor2.p, inductor2.n) annotation(
        Line(points = {{0, -90}, {0, -86}, {-7, -86}}, color = {0, 0, 255}));
      connect(capacitor2.n, resistor14.p) annotation(
        Line(points = {{0, -110}, {11, -110}}, color = {0, 0, 255}));
      connect(resistor13.n, inductor8.p) annotation(
        Line(points = {{35, -86}, {43, -86}}, color = {0, 0, 255}));
      connect(and2.y, switch2.control) annotation(
        Line(points = {{315, -184}, {324, -184}, {324, -121}, {399, -121}, {399, 10}, {413, 10}}, color = {255, 0, 255}));
      connect(and1.y, idealClosingSwitch.control) annotation(
        Line(points = {{389, -184}, {405, -184}, {405, -48}, {416, -48}}, color = {255, 0, 255}));
      connect(and3.y, idealClosingSwitch1.control) annotation(
        Line(points = {{437, -201}, {449, -201}, {449, -141}, {412, -141}, {412, -106}, {416, -106}}, color = {255, 0, 255}));
      connect(realExpression6.y, dsogi_pll.V_abc[1]) annotation(
        Line(points = {{-144, 175}, {-135, 175}, {-135, 160}, {-118, 160}}, color = {0, 0, 127}));
      connect(realExpression4.y, dsogi_pll.V_abc[2]) annotation(
        Line(points = {{-144, 159}, {-131, 159}, {-131, 160}, {-118, 160}}, color = {0, 0, 127}));
      connect(realExpression5.y, dsogi_pll.V_abc[3]) annotation(
        Line(points = {{-144, 143}, {-135, 143}, {-135, 160}, {-118, 160}}, color = {0, 0, 127}));
      connect(dsogi_pll.d, vsd.u1) annotation(
        Line(points = {{-96, 169}, {-82, 169}, {-82, 145}, {-49, 145}}, color = {0, 0, 127}));
      connect(dsogi_pll.q, vsq.u1) annotation(
        Line(points = {{-96, 163}, {-63, 163}, {-63, 121}, {-50, 121}}, color = {0, 0, 127}));
      connect(realExpression7.y, innerControl_PIR.vd) annotation(
        Line(points = {{386.7, 192.5}, {386.7, 193.5}, {408, 193.5}, {408, 183}}, color = {0, 0, 127}));
      connect(realExpression10.y, innerControl_PIR.Id) annotation(
        Line(points = {{385.9, 164}, {394.4, 164}, {394.4, 171}, {408, 171}}, color = {0, 0, 127}));
      connect(realExpression11.y, innerControl_PIR.Iq) annotation(
        Line(points = {{386.95, 145}, {397.45, 145}, {397.45, 165}, {408, 165}}, color = {0, 0, 127}));
      connect(innerControl_PIR.vdref, division.u1) annotation(
        Line(points = {{443, 178}, {445.05, 178}, {445.05, 185.18}, {459.55, 185.18}}, color = {0, 0, 127}));
      connect(innerControl_PIR.vqref, division1.u1) annotation(
        Line(points = {{443, 162}, {443, 148.37}, {457.55, 148.37}}, color = {0, 0, 127}));
      connect(innerControl_PIR.vq, realExpression8.y) annotation(
        Line(points = {{408, 156}, {404.45, 156}, {404.45, 117.55}, {353.45, 117.55}}, color = {0, 0, 127}));
      connect(realExpression3.y, Correntes.theta) annotation(
        Line(points = {{258.05, 137}, {267, 137}, {267, 149}}, color = {0, 0, 127}));
      connect(vsa.p, inductor6.n) annotation(
        Line(points = {{69, 52}, {60, 52}, {60, 33.5}}, color = {0, 0, 255}));
      connect(vsb.p, inductor7.n) annotation(
        Line(points = {{65, -5}, {60, -5}, {60, -25}}, color = {0, 0, 255}));
      connect(vscv.p, inductor8.n) annotation(
        Line(points = {{69, -61}, {63, -61}, {63, -86}}, color = {0, 0, 255}));
      connect(vsc.Vta, resistor.p) annotation(
        Line(points = {{-123.7, 34.8}, {-91.7, 34.8}, {-91.7, 33.8}, {-59.7, 33.8}}, color = {0, 0, 255}));
      connect(vsc.Vtb, resistor1.p) annotation(
        Line(points = {{-123.7, 27}, {-82.7, 27}, {-82.7, -26}, {-56.7, -26}}, color = {0, 0, 255}));
      connect(vsc.Vtc, resistor2.p) annotation(
        Line(points = {{-123.7, 19.2}, {-92.7, 19.2}, {-92.7, -85.8}, {-54.7, -85.8}}, color = {0, 0, 255}));
      connect(inductor6.n, ia.p) annotation(
        Line(points = {{60, 33.5}, {93, 33.5}, {93, 33}}, color = {0, 0, 255}));
      connect(inductor7.n, ib.p) annotation(
        Line(points = {{60, -25}, {95, -25}, {95, -27}}, color = {0, 0, 255}));
      connect(inductor8.n, ic.p) annotation(
        Line(points = {{63, -86}, {79.5, -86}, {79.5, -87}, {97, -87}}, color = {0, 0, 255}));
      connect(realExpression9.y, innerControl_PIR.Id_ref) annotation(
        Line(points = {{387.8, 177.5}, {408, 177.5}, {408, 177}}, color = {0, 0, 127}));
      connect(realExpression12.y, innerControl_PIR.Iqref) annotation(
        Line(points = {{387.6, 133.5}, {400.6, 133.5}, {400.6, 160}, {408, 160}}, color = {0, 0, 127}));
      connect(pca.p, ipca.p) annotation(
        Line(points = {{368, 14}, {345, 14}, {345, 38}}, color = {0, 0, 255}));
      connect(pcb.p, ipcb.p) annotation(
        Line(points = {{371, -45}, {347, -45}, {347, -22}}, color = {0, 0, 255}));
      connect(pcc.p, ipcc.p) annotation(
        Line(points = {{366, -107}, {349, -107}, {349, -82}}, color = {0, 0, 255}));
      connect(realExpression47.y, ipcb_pu.u1) annotation(
        Line(points = {{645, -161}, {683, -161}, {683, -159}}, color = {0, 0, 127}));
      connect(realExpression49.y, ipca_pu.u1) annotation(
        Line(points = {{645, -145}, {660, -145}, {660, -131}, {683, -131}}, color = {0, 0, 127}));
      connect(constant17.y, ipcc_pu.u2) annotation(
        Line(points = {{665, -208}, {673, -208}, {673, -193}, {683, -193}}, color = {0, 0, 127}));
      connect(constant17.y, ipcb_pu.u2) annotation(
        Line(points = {{665, -208}, {673, -208}, {673, -167}, {683, -167}}, color = {0, 0, 127}));
      connect(realExpression48.y, ipcc_pu.u1) annotation(
        Line(points = {{645, -177}, {676, -177}, {676, -185}, {683, -185}}, color = {0, 0, 127}));
      connect(constant17.y, ipca_pu.u2) annotation(
        Line(points = {{665, -208}, {673, -208}, {673, -139}, {683, -139}}, color = {0, 0, 127}));
      connect(Pramp.y, Controle_potencia.P_ref) annotation(
        Line(points = {{33, 206}, {69, 206}, {69, 146}, {77, 146}}, color = {0, 0, 127}));
      connect(Pramp.y, Pref_pu.u1) annotation(
        Line(points = {{33, 206}, {103, 206}, {103, 196}, {110, 196}}, color = {0, 0, 127}));
      connect(constant18.y, Qref_pu.u2) annotation(
        Line(points = {{104, 90}, {118, 90}, {118, 105}, {127, 105}}, color = {0, 0, 127}));
      connect(constant12.y, pca_pu.u2) annotation(
        Line(points = {{862, -151}, {867, -151}, {867, -142}, {876, -142}}, color = {0, 0, 127}));
      connect(realExpression20.y, pca_pu.u1) annotation(
        Line(points = {{865, -131}, {870, -131}, {870, -134}, {876, -134}}, color = {0, 0, 127}));
      connect(Qref_pu.u1, Controle_potencia.Q_ref) annotation(
        Line(points = {{126.6, 113.2}, {70.6, 113.2}, {70.6, 127.2}, {76.6, 127.2}}, color = {0, 0, 127}));
      connect(constantVoltage.n, constantVoltage1.p) annotation(
        Line(points = {{-223, 8}, {-223, -6}}, color = {0, 0, 255}));
      connect(constantVoltage1.n, vsc.pin_n) annotation(
        Line(points = {{-223, -26}, {-197, -26}, {-197, 20.5}, {-152, 20.5}}, color = {0, 0, 255}));
      connect(currentSensor.n, constantVoltage.n) annotation(
        Line(points = {{-230, 1}, {-223, 1}, {-223, 8}}, color = {0, 0, 255}));
      connect(ground.p, currentSensor.p) annotation(
        Line(points = {{-248, 0}, {-238, 0}, {-238, 1}}, color = {0, 0, 255}));
      connect(resistor9.n, filtro_a.p) annotation(
        Line(points = {{27, 10}, {33, 10}}, color = {0, 0, 255}));
      connect(resistor14.n, filtro_c.p) annotation(
        Line(points = {{31, -110}, {36, -110}}, color = {0, 0, 255}));
      connect(resistor12.n, filtro_b.p) annotation(
        Line(points = {{28, -50}, {33, -50}}, color = {0, 0, 255}));
      connect(RG2.n, LG2.p) annotation(
        Line(points = {{399, -22}, {406, -22}, {406, -21}, {412, -21}}, color = {0, 0, 255}));
      connect(sineVoltage1.p, LG2.n) annotation(
        Line(points = {{445, -21}, {432, -21}}, color = {0, 0, 255}));
      connect(LG3.n, sineVoltage2.p) annotation(
        Line(points = {{433, -80}, {445, -80}}, color = {0, 0, 255}));
      connect(LG3.p, RG3.n) annotation(
        Line(points = {{413, -80}, {404, -80}, {404, -82}, {399, -82}}, color = {0, 0, 255}));
      connect(d_D_Transformer.pin3, ia.n) annotation(
        Line(points = {{270, -13}, {257, -13}, {257, 33}, {105, 33}}, color = {0, 0, 255}));
      connect(d_D_Transformer.pin4, ib.n) annotation(
        Line(points = {{270, -27}, {107, -27}}, color = {0, 0, 255}));
      connect(d_D_Transformer.pin5, ic.n) annotation(
        Line(points = {{271, -40}, {252, -40}, {252, -87}, {109, -87}}, color = {0, 0, 255}));
      connect(d_D_Transformer.pin, ipca.p) annotation(
        Line(points = {{328, -14}, {333, -14}, {333, 38}, {345, 38}}, color = {0, 0, 255}));
      connect(d_D_Transformer.pin1, ipcb.p) annotation(
        Line(points = {{328, -27}, {339, -27}, {339, -22}, {347, -22}}, color = {0, 0, 255}));
      connect(d_D_Transformer.pin2, ipcc.p) annotation(
        Line(points = {{327, -40}, {338, -40}, {338, -82}, {349, -82}}, color = {0, 0, 255}));
      connect(filtro_a.n, filtro_b.n) annotation(
        Line(points = {{45, 10}, {45, -50}}, color = {0, 0, 255}));
      connect(filtro_b.n, filtro_c.n) annotation(
        Line(points = {{45, -50}, {48, -50}, {48, -110}}, color = {0, 0, 255}));
      connect(filtro_c.n, filtro_terra.p) annotation(
        Line(points = {{48, -110}, {57, -110}, {57, -111}}, color = {0, 0, 255}));
      connect(filtro_terra.n, ground14.p) annotation(
        Line(points = {{69, -111}, {76, -111}}, color = {0, 0, 255}));
      connect(Qramp.y, Controle_potencia.Q_ref) annotation(
        Line(points = {{43, 102}, {64, 102}, {64, 127}, {77, 127}}, color = {0, 0, 127}));
      connect(realExpression50.y, vsa_pu.u1) annotation(
        Line(points = {{935, -131}, {940, -131}, {940, -134}, {946, -134}}, color = {0, 0, 127}));
      connect(constant19.y, vsa_pu.u2) annotation(
        Line(points = {{932, -151}, {937, -151}, {937, -142}, {946, -142}}, color = {0, 0, 127}));
      connect(vsab.p, ia.n) annotation(
        Line(points = {{115, 14}, {114, 14}, {114, 33}, {105, 33}}, color = {0, 0, 255}));
      connect(vsab.n, ib.n) annotation(
        Line(points = {{115, 0}, {116, 0}, {116, -27}, {107, -27}}, color = {0, 0, 255}));
      connect(vsbc.p, ib.n) annotation(
        Line(points = {{116, -52}, {117, -52}, {117, -27}, {107, -27}}, color = {0, 0, 255}));
      connect(vsbc.n, ic.n) annotation(
        Line(points = {{116, -66}, {115, -66}, {115, -87}, {109, -87}}, color = {0, 0, 255}));
      connect(vsac.p, ia.n) annotation(
        Line(points = {{123, -9}, {123, 33}, {105, 33}}, color = {0, 0, 255}));
      connect(vsac.n, ic.n) annotation(
        Line(points = {{123, -23}, {123, -87}, {109, -87}}, color = {0, 0, 255}));
      connect(resistor19.n, Loadgb.p) annotation(
        Line(points = {{223, -74}, {200, -74}, {200, -73}}, color = {0, 0, 255}));
      connect(Rloadb.p, d_D_Transformer.pin4) annotation(
        Line(points = {{200, -34}, {200, -27}, {270, -27}}, color = {0, 0, 255}));
      connect(Rloadc.p, d_D_Transformer.pin5) annotation(
        Line(points = {{145, -31}, {183, -31}, {183, -87}, {252, -87}, {252, -40}, {271, -40}}, color = {0, 0, 255}));
      connect(Loadgb.n, d_D_Transformer.pin5) annotation(
        Line(points = {{200, -85}, {242, -85}, {242, -40}, {271, -40}}, color = {0, 0, 255}));
      connect(Loadgc.n, d_D_Transformer.pin3) annotation(
        Line(points = {{145, -83}, {133, -83}, {133, -13}, {270, -13}}, color = {0, 0, 255}));
      connect(constant14.y, vsb_pu.u2) annotation(
        Line(points = {{989, -92}, {994, -92}, {994, -83}, {1003, -83}}, color = {0, 0, 127}));
      connect(realExpression31.y, pcb_pu.u1) annotation(
        Line(points = {{922, -72}, {927, -72}, {927, -75}, {933, -75}}, color = {0, 0, 127}));
      connect(realExpression30.y, vsb_pu.u1) annotation(
        Line(points = {{992, -72}, {997, -72}, {997, -75}, {1003, -75}}, color = {0, 0, 127}));
      connect(constant3.y, pcb_pu.u2) annotation(
        Line(points = {{919, -92}, {924, -92}, {924, -83}, {933, -83}}, color = {0, 0, 127}));
      connect(realExpression38.y, vscv_pu.u1) annotation(
        Line(points = {{837, -260}, {842, -260}, {842, -263}, {848, -263}}, color = {0, 0, 127}));
      connect(realExpression32.y, pcc_pu.u1) annotation(
        Line(points = {{767, -260}, {772, -260}, {772, -263}, {778, -263}}, color = {0, 0, 127}));
      connect(constant20.y, pcc_pu.u2) annotation(
        Line(points = {{764, -280}, {769, -280}, {769, -271}, {778, -271}}, color = {0, 0, 127}));
      connect(constant15.y, vscv_pu.u2) annotation(
        Line(points = {{834, -280}, {839, -280}, {839, -271}, {848, -271}}, color = {0, 0, 127}));
      connect(currentSensor2.n, d_D_Transformer.pin4) annotation(
        Line(points = {{196, -21}, {259, -21}, {259, -27}, {270, -27}}, color = {0, 0, 255}));
      connect(currentSensor1.p, d_D_Transformer.pin3) annotation(
        Line(points = {{196, 31}, {252, 31}, {252, -13}, {270, -13}}, color = {0, 0, 255}));
      connect(inductor5.n, resistor20.p) annotation(
        Line(points = {{166, -58}, {166, -62}}, color = {0, 0, 255}));
      connect(inductor3.n, resistor15.p) annotation(
        Line(points = {{217, 5}, {217, 2}}, color = {0, 0, 255}));
      connect(inductor4.n, resistor19.p) annotation(
        Line(points = {{223, -50}, {223, -54}}, color = {0, 0, 255}));
      connect(inductor3.p, currentSensor1.n) annotation(
        Line(points = {{217, 25}, {204, 25}, {204, 15}, {196, 15}, {196, 19}}, color = {0, 0, 255}));
      connect(resistor15.n, currentSensor2.p) annotation(
        Line(points = {{217, -18}, {202, -18}, {202, -2}, {196, -2}, {196, -9}}, color = {0, 0, 255}));
      connect(inductor5.p, Rloadc.n) annotation(
        Line(points = {{166, -38}, {154, -38}, {154, -47}, {145, -47}, {145, -43}}, color = {0, 0, 255}));
      connect(resistor20.n, Loadgc.p) annotation(
        Line(points = {{166, -82}, {166, -85}, {151, -85}, {151, -63}, {145, -63}, {145, -71}}, color = {0, 0, 255}));
      connect(inductor4.p, Rloadb.n) annotation(
        Line(points = {{223, -30}, {210, -30}, {210, -46}, {200, -46}}, color = {0, 0, 255}));
      annotation(
        Diagram(coordinateSystem(extent = {{-150, 330}, {790, -240}}, grid = {1, 1})),
        Icon(coordinateSystem(extent = {{-120, -120}, {480, 160}}, grid = {1, 1})));
    end SAEB1_S2_2;
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
