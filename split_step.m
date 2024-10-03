solverversion = "GPE solver v1.4, RZ 2019";
Print[solverversion];
Off[General::munfl ];
algorithm = "firstOrder";
evolution[type_String: "ite", steps_: 1, monitorSteps_Integer: 1, 
   plotSteps_Integer: 0] := Module[{s, calcSteps},
   Assert[Mod[plotSteps, monitorSteps] == 0];
   monitor;
   If[plotSteps > 0, plot];
   s = 0;
   While[s < steps,
    calcSteps = Min[{steps - s, monitorSteps}];
    evolutionCmd[type, calcSteps, algorithm];
    s = s + calcSteps;
    monitor;
    If[plotSteps > 0 && (Mod[s, plotSteps] == 0 || s == steps), 
     plot];
    ];
   ];
evolve[type_String: "ite", tmax_, nrmonitor_: 1, nrplot_: 0] := 
  Module[{steps = tmax/timeStep},
   evolution[type, Round[steps], Round[steps/nrmonitor], 
    If[nrplot > 0, Round[steps/nrplot], 0]]
   ];
defaultscale = 5;
HarmonicTrap[Omega_List] := Module[{freq0, time0, length0},
   omega = Omega;
   (* rescaling factors *)
   freq0 = Omega[[-1]];
   time0 = 1/freq0;
   length0 = Sqrt[1/freq0];
   ho = Sqrt[1/omega]/length0;
   max = defaultscale ho;
   min = -max;
   cho = omega^2/freq0^2;
   hofnc[r_List] := cho.r^2/2;
   potfnc = hofnc;
   reset;
   ];
HarmonicTrap3D[OmegaX_: 1, OmegaY_: 1, OmegaZ_: 1] := 
  HarmonicTrap[{OmegaX, OmegaY, OmegaZ}];
HarmonicTrap2D[OmegaX_: 1, OmegaY_: 1] := 
  HarmonicTrap[{OmegaX, OmegaY}];
HarmonicTrap1D[OmegaX_: 1] := HarmonicTrap[{OmegaX}];
NoPotential[scale_: defaultscale] := Module[{},
   max = Table[scale, {dim}];
   min = -max;
   ho = Table[1, {dim}];
   potfnc = 0 &;
   reset;
   ];
gridReport := Grid[{
    {"sizes", "dim", "len", "max", "min", "delta", "deltak", "maxk"},
    {sizes, dim, len, max, min, delta, deltak, maxk}}];
dim := Length[sizes];
len := Times @@ sizes;
indexes[dir_] := Table[{n}, {n, sizes[[dir]]}];
mesh[dir_] := Table[N[xx[dir, n]], {n, sizes[[dir]]}];
kmesh[dir_] := Table[N[kk[dir, n]], {n, sizes[[dir]]}];
indexes[] := 
  Outer[List, Sequence @@ Table[indexes[dir], {dir, dim}]];
mesh[] := Outer[List, Sequence @@ Table[mesh[dir], {dir, dim}]];
kmesh[] := Outer[List, Sequence @@ Table[kmesh[dir], {dir, dim}]];
meshX := mesh[1];
meshY := mesh[2];
meshZ := mesh[3];
xx[n_, i_] := min[[n]] + (i - 1) delta[[n]];
kk[n_, i_] /; (i - 1) < sizes[[n]]/2 := (i - 1) deltak[[n]];
kk[n_, i_] /; (i - 1) >= 
    sizes[[n]]/2 := -maxk[[n]] + ((i - 1) - sizes[[n]]/2)*
    deltak[[n]];
dV := Times @@ delta;
Gaussian[s_List] := Module[{factor, fnc},
   factor = Pi^(-dim/4) (Times @@ s)^(-1/2);
   fnc[r_List] := N[factor Exp[-Total[r^2/(2 s^2)]]];
   Map[fnc, mesh[], {dim}]
   ];
initPK := potentialKinetic = Map[N[#.#/2] &, kmesh[], {dim}];
initPE := potentialExternal = Map[N[potfnc[#]] &, mesh[], {dim}];
initialize[initwf_: True] := Module[{},
   delta = (max - min)/(sizes - 1);
   deltak = 2 Pi/(sizes delta);
   maxk = Pi/delta;
   If[initwf,
    wf = Gaussian[ho];
    wf /= Sqrt[norm];
    ];
   initPK;
   initPE;
   ];
newPotential[pot_] := {potfnc = pot; initPE};
free := 0 &;
evolcoef["rte"] := -I;
evolcoef["ite"] := -1;
normalization["ite"] := Sqrt[norm];
normalization["rte"] := 1;
threeBodyLosses = False;
ft[wf_] := Fourier[wf, FourierParameters -> {1, -1}];
ift[wf_] := InverseFourier[wf, FourierParameters -> {1, -1}];
evolutionCmd[type_String, steps_Integer, "fourthOrder"] := 
  Module[{dt, evolutionCoefficient, matAc1, matAc2, matBd1, matBd2, 
    density, step,
    x0 = -2^(1/3)/(2 - 2^(1/3)), x1 = 1/(2 - 2^(1/3)), c1, c2, d1, d2},
   c1 = x1/2; c2 = (x0 + x1)/2; d1 = x1; d2 = x0;
   dt = timeStep;
   evolutionCoefficient = evolcoef[type] dt;
   matAc1 = Exp[potentialKinetic evolutionCoefficient c1];
   matAc2 = Exp[potentialKinetic evolutionCoefficient c2];
   matBd1 = Exp[potentialExternal evolutionCoefficient d1];
   matBd2 = Exp[potentialExternal evolutionCoefficient d2];
   Do[
    wf = ft[wf];
    wf *= matAc1;
    wf = ift[wf];
    
    density = Abs[wf]^2;
    wf *= 
     Exp[density contactInteractionFactor evolutionCoefficient d1];
    wf *= matBd1;
    
    wf = ft[wf];
    wf *= matAc2;
    wf = ift[wf];
    
    density = Abs[wf]^2;
    wf *= 
     Exp[density contactInteractionFactor evolutionCoefficient d2];
    wf *= matBd2;
    
    wf = ft[wf];
    wf *= matAc2; (* c3 *)
    wf = ift[wf];
    
    density = Abs[wf]^2;
    wf *= 
     Exp[density contactInteractionFactor evolutionCoefficient d1]; (* 
    d3 *)
    wf *= matBd1; (* d3 *)
    
    wf = ft[wf];
    wf *= matAc1; (* c4 *)
    wf = ift[wf];
    
    wf /= normalization[type];
    time += timeStep,
    {step, steps}];
   totalSteps += steps;
   steps
   ];
evolutionCmd[type_String, steps_Integer, "secondOrder"] := 
  Module[{dt, evolutionCoefficient, matAhalf, density, step},
   dt = timeStep;
   evolutionCoefficient = evolcoef[type] dt;
   matAhalf = Exp[potentialKinetic evolutionCoefficient/2];
   Do[
    wf = ft[wf];
    wf *= matAhalf;
    wf = ift[wf];
    
    density = Abs[wf]^2;
    wf *= 
     Exp[evolutionCoefficient (density contactInteractionFactor + 
         potentialExternal )];
    
    wf = ft[wf];
    wf *= matAhalf;
    wf = ift[wf];
    
    wf /= normalization[type];
    time += timeStep,
    {step, steps}];
   totalSteps += steps;
   steps
   ];
(* Three-body losses only implemented for first-order approach. *)

evolutionCmd[type_String, steps_Integer, "firstOrder"] := 
  Module[{dt, evolutionCoefficient, matA, matB, density, step},
   dt = timeStep;
   evolutionCoefficient = evolcoef[type] dt;
   matA = Exp[potentialKinetic evolutionCoefficient];
   matB = Exp[potentialExternal evolutionCoefficient];
   Do[
    density = Abs[wf]^2;
    wf *= Exp[density contactInteractionFactor evolutionCoefficient];
    If[threeBodyLosses === True,
     wf *= Exp[density^2 threeBodyLossesFactor evolutionCoefficient];
     ];
    wf *= matB;
    
    wf = ft[wf];
    wf *= matA;
    wf = ift[wf];
    
    wf /= normalization[type];
    time += timeStep,
    {step, steps}];
   totalSteps += steps;
   steps
   ];
intg[wf_, f_] := Chop[Total[Conjugate[wf] f wf, dim] dV];
expv[fnc_] := Module[{f},
   f = Map[fnc, mesh[], {dim}];
   intg[wf, f]
   ];
norm := Chop[Total[Abs[wf]^2, dim] dV];
mom[n_, k_] := expv[#[[n]]^k &];
mean[n_] := mom[n, 1];
meanX := mean[1];
meanY := mean[2];
meanZ := mean[3];
sigma[n_] := Sqrt[mom[n, 2]];
sigmaX := sigma[1];
sigmaY := sigma[2];
sigmaZ := sigma[3];
dir[1] := "X"; dir[2] := "Y"; dir[3] := "Z";
moments := Join[Array[mean, dim], Array[sigma, dim], {norm}];
momentslegend := 
  Join[Table["mean" <> dir[i], {i, dim}], 
   Table["sigma" <> dir[i], {i, dim}], {"norm"}];
energyKinetic := Module[{wfk = ft[wf]},
   intg[wfk, potentialKinetic]/len
   ];
energyPotential := intg[wf, potentialExternal];
energyContact := Module[{density = Abs[wf]^2},
   contactInteractionFactor intg[wf, density]/2
   ];
energyInt := energyContact;
energy := energyKinetic + energyPotential + energyInt;
pow = 2; (* h.o. *)

virial := 2 energyKinetic - pow energyPotential + dim energyInt;
energies := {energy, energyKinetic, energyPotential, energyContact, 
   virial};
energieslegend := {"total", "kinetic", "potential", "contact", 
   "virial"};
misc := {totalSteps, Abs[ psi[origin] ]};
misclegend := {"steps", "A0"};
myTotal[l_, d : {1, 2}] := Total[l, d];
myTotal[l_, d : {1, 3}] := Total[Total[l, {1}], {2}];
myTotal[l_, d : {2, 3}] := Total[l, d];
myTotal[l_, {n_}] := Total[l, {n}];
myTotal[l_, {}] := l;
projection[d_] := Module[{compl, vol},
   compl = Complement[Range[dim], {d}];
   vol = Times @@ Part[delta, compl];
   Transpose[{mesh[d], myTotal[Abs[wf]^2, compl] vol}]
   ];
projectionX := projection[1];
projectionY := projection[2];
projectionZ := projection[3];
projections := Table[projection[i], {i, dim}];
plotProjectionsOpts[opts___] := 
  ListPlot[Array[projection, {dim}], Joined  True, PlotRange  All, 
   Sequence[opts]];
plotProjections := plotProjectionsOpts[];
reset := {
   time = 0;
   totalSteps = 0;
   lenergies = {};
   lmoments = {};
   lmisc = {};
   lall = {};
   lprojections = {};
   lplots = {};
   };
calcEnergies := Join[{time}, energies];
calcMoments := Join[{time}, moments];
calcMisc := Join[{time}, misc];
calcAll := {{time}, Transpose[{energieslegend, energies}], 
   Transpose[{momentslegend, moments}], Transpose[{misclegend, misc}]};
calcProjections := Join[{time}, projections];
calcAllTable := 
  Grid[{{"time", "energies", "moments", "misc"}, 
    Map[MatrixForm, calcAll]}];
myPrint = PrintTemporary;
monitor := Module[{l1, l2, l3},
   myPrint[Grid[{Map[Grid[{#}, Alignment  Left] &, calcAll]}]];
   AppendTo[lenergies, l1 = calcEnergies];
   AppendTo[lmoments, l2 = calcMoments];
   AppendTo[lmisc, l3 = calcMisc];
   AppendTo[lall, Join[{time}, Rest[l1], Rest[l2], Rest[l3]]];
   AppendTo[lprojections, calcProjections];
   ];
interp[wf_] := Module[{tab},
   tab = MapIndexed[{#1, Part[wf, Sequence @@ #2]} &, mesh[], {dim}];
   tab = Flatten[tab, dim - 1];
   Interpolation[tab]
   ];
psi := interp[wf];
origin := Sequence @@ Table[0, {dim}];
plot := AppendTo[lplots, {time, plotProjections, psi}];
reportEnergies := 
  TableForm[Prepend[lenergies, Prepend[energieslegend, "time"]]];
reportMoments := 
  TableForm[Prepend[lmoments, Prepend[momentslegend, "time"]]];
reportMisc := TableForm[Prepend[lmisc, Prepend[misclegend, "time"]]];
legend := Join[energieslegend, momentslegend, misclegend];
report := TableForm[Prepend[lall, Prepend[legend, "time"]]];
getlist[l_, legend_, key_] := 
  Module[{pos = FirstPosition[legend, key][[1]]},
   l[[All, {1, pos + 1}]]
   ];
listM[key_] /; MemberQ[legend, key] := getlist[lall, legend, key];
plotM[key_] /; MemberQ[legend, key] := 
  ListPlot[listM[key], Joined -> True, PlotRange -> All, 
   AxesLabel -> {"t", key}];