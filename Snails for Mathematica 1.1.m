(* ::Package:: *)

(* ::Title:: *)
(*Snail functions for Mathematica*)


(* ::Text:: *)
(*Version 1.1*)
(*15 October 2025*)
(**)
(*This notebook contains functions for simulating and processing virtual gastropod (and other mollusc) shells.  *)
(**)
(*P. David Polly*)
(*pdpolly@pollylab.org*)
(*https://pollylab.org/*)


(*  This function prints the version number for the installed verison of this package.  *)
VersionNumber="1.1";
VersionDate="15 October 2025";
SnailsVersion[]:=Print["Snails for Mathematica "<>VersionNumber<>"\n(c) P. David Polly, "<>VersionDate<>"\n"];
SnailsVersion[]


(* This function generates a three-dimensional rendered shell using the functions 
	developed by David Raup (1966). 

	Usage:  RaupCoil3D[shape, W, T, D, TurnNum]
			where shape is a matrix of 3D semilandmarks that circumscribe the aperture
			W is Raup's whorl expansion rate parameter
			T is Raup's vertical translation rate parameter
			D is Raup's parameter for distance of aperture from columella
			TurnNum is the desired number of whorls 
			
	Created 16 March 2016. 
	Updated 5 February 2023 with adjustment to aperture position relative to coiling. 
*)

RaupCoil3D[S_,W_,T_,Di_, turns_,CoordsOnly_:False]:=Quiet[Block[{r,z,x,CoiledCoords,\[Theta],Sprime,S3D,yaspect,rc},

Sprime=#-Mean[S]&/@S;
Sprime=(Sprime/(Max[Sprime[[1;;,1]]]-Min[Sprime[[1;;,1]]]));
Sprime[[1;;,1]]=Sprime[[1;;,1]]-Min[Sprime[[1;;,1]]]+Di;
yaspect=Max[Sprime[[1;;,3]]]-Min[Sprime[[1;;,3]]];
rc=Mean[Sprime[[1;;,1]]];
S3D=Sprime;
CoiledCoords=Table[Table[
r=S3D[[x,1]]*(W^(\[Theta]/(2Pi)));
z=(S3D[[x,3]]*(W^(\[Theta]/(2Pi))))+(T*(W^(\[Theta]/(2Pi))-0.5));
{r*Cos[\[Theta]],r*Sin[\[Theta]]+S3D[[x,2]]*(W^(\[Theta]/(2Pi))),z}
,{x,Length[S3D]}]
,{\[Theta],0,turns*2Pi,10Degree}];
If[CoordsOnly==True,Return[CoiledCoords],
Return[
Graphics3D[
Table[{
Table[{RGBColor[0.997681, 0.980789, 0.842313],EdgeForm[None],Opacity[1],Polygon[{CoiledCoords[[x,y]],CoiledCoords[[x-1,y]],CoiledCoords[[x-1,y+1]],CoiledCoords[[x,y+1]]}]},{y,Length[CoiledCoords[[x]]]-1}],
{GrayLevel[0.5],Opacity[0],AbsolutePointSize[1.75],Point[CoiledCoords[[x]]]}},
{x,2,Length[CoiledCoords]}],Boxed->False]]];
]];



(* This function generates a three-dimensional rendered shell using the functions 
	developed by David Raup (1966) as a three-dimensional mesh suitable for printing
	with a 3D preinter.  It is nearly idential ato RaupCoil3D execpt it creates a shell
	with two layers that is capped at apex and aperture. 

	Usage:  RaupCoil3DForPrint[shape, W, T, D, TurnNum]
			where shape is a matrix of 3D semilandmarks that circumscribe the aperture
			W is Raup's whorl expansion rate parameter
			T is Raup's vertical translation rate parameter
			D is Raup's parameter for distance of aperture from columella
			TurnNum is the desired number of whorls 
			
	Created 16 March 2016.  
*)

RaupCoil3DForPrint[S_,W_,T_,Di_, turns_]:=Quiet[Block[{r,z,x,CoiledCoords,\[Theta],Sprime,S3D,yaspect,rc},

Sprime=#-Mean[S]&/@S;
Sprime=(Sprime/(Max[Sprime[[1;;,1]]]-Min[Sprime[[1;;,1]]]));
Sprime[[1;;,1]]=Sprime[[1;;,1]]-Min[Sprime[[1;;,1]]]+Di;
yaspect=Max[Sprime[[1;;,3]]]-Min[Sprime[[1;;,3]]];
rc=Mean[Sprime[[1;;,1]]];
S3D=Sprime;
CoiledCoords=Table[Table[
r=S3D[[x,1]]*(W^(\[Theta]/(2Pi)));
z=(S3D[[x,3]]*(W^(\[Theta]/(2Pi))))+rc*T*yaspect*((W^(\[Theta]/(2Pi)))-1);
{r*Cos[\[Theta]],r*Sin[\[Theta]]+S3D[[x,2]]*(W^(\[Theta]/(2Pi))),z}
,{x,Length[S3D]}]
,{\[Theta],0,turns*2Pi,10Degree}];
Return[CoiledCoords];
]];


(* This function tilts the aperture shape at a specified number of degrees. 

	Usage:  TiltAperture[shape, angle]
			where shape is a matrix of 3D semilandmarks that circumscribe the aperture
			and angle is an angle from vertical given in degrees
			
	Created 16 March 2016.  
*)

TiltAperture[aperture_,angle_]:=Module[{tiltedaperture},
tiltedaperture=aperture . RotationMatrix[angle Degree,{1,0,0}];
tiltedaperture[[1;;,2]]=tiltedaperture[[1;;,2]]+Mean[tiltedaperture[[1;;,2]]];
Return[tiltedaperture];
];


(* This function generates three-dimensional elliptical Fourier coefficiences from 
	a closed curve (like an aperture shape).  It returns Fourier cofficients for a
	single shape.

	Usage:  EllipticalFourierCoefficients[points]
			where points is a matrix of 3D semilandmarks that circumscribe the aperture
			
	Created 16 March 2016.  
*)

EllipticalFourierCoefficients[standardizedoutline_]:=Module[{deltaXY,deltaT,TotalT,FourierCoeffs},
deltaXY=Table[standardizedoutline[[k+1]]-standardizedoutline[[k]],{k,Length[standardizedoutline]-1}];
deltaT=Sqrt[Plus@@(#^2)]&/@deltaXY;
TotalT=Plus@@deltaT;
FourierCoeffs=Table[{TotalT/(2*(n^2)*(Pi^2))*(Plus@@Table[deltaXY[[k,1]]/deltaT[[k]]*(Cos[(2Pi*n*Plus@@deltaT[[1;;k+1]])/TotalT]-Cos[(2Pi*n*Plus@@deltaT[[1;;k]])/TotalT]),{k,Length[deltaXY]-1}]),TotalT/(2*(n^2)*(Pi^2))*(Plus@@Table[deltaXY[[k,1]]/deltaT[[k]]*(Sin[(2Pi*n*Plus@@deltaT[[1;;k+1]])/TotalT]-Sin[(2Pi*n*Plus@@deltaT[[1;;k]])/TotalT]),{k,Length[deltaXY]-1}]),TotalT/(2*(n^2)*(Pi^2))*(Plus@@Table[deltaXY[[k,2]]/deltaT[[k]]*(Cos[(2Pi*n*Plus@@deltaT[[1;;k+1]])/TotalT]-Cos[(2Pi*n*Plus@@deltaT[[1;;k]])/TotalT]),{k,Length[deltaXY]-1}]),TotalT/(2*(n^2)*(Pi^2))*(Plus@@Table[deltaXY[[k,2]]/deltaT[[k]]*(Sin[(2Pi*n*Plus@@deltaT[[1;;k+1]])/TotalT]-Sin[(2Pi*n*Plus@@deltaT[[1;;k]])/TotalT]),{k,Length[deltaXY]-1}])},{n,Length[deltaXY]/2}];
Return[{FourierCoeffs,deltaT,TotalT}];
];


(* This function generates a shape in  X, Y coordinates from elliptical Fourier coefficients.

	Usage:  HarmonicReconstruction[coefficients, deltaT, T, n]
			where coefficients is a matrix of Fourier coefficients
			deltaT is the offset 
			T is the number
			n is 
			
	Created 16 March 2016.  
*)
HarmonicReconstruction[coeffs_,deltaT_,T_,n_]:=Module[{xy},
xy=Table[Flatten[{Plus@@Table[coeffs[[j,1]]*Cos[(2Pi*j*Plus@@deltaT[[1;;i]])/T]+coeffs[[j,2]]*Sin[(2Pi*j*Plus@@deltaT[[1;;i]])/T],{j,n}],Plus@@Table[coeffs[[j,3]]*Cos[(2Pi*j*Plus@@deltaT[[1;;i]])/T]+coeffs[[j,4]]*Sin[(2Pi*j*Plus@@deltaT[[1;;i]])/T],{j,n}]}],{i,Length[deltaT]}];
Return[xy]
];


(* This function saves a 3D shell suitable for printing as an OBJ mesh file.

	Usage:  SaveShell[shellobject, parameters, title, path]
			where shellobject is a shell generated with the RaupCoil3DForPrint function
			parameters are the shell coiling parameters used to generate the shell(which are included as metadata in the OBJ file)
			title is a title for the shell object (included as metadata in the OBJ file)
			path is the path to the location where the file will be saved.
			
	Created 16 March 2016.  
*)

SaveShell[myShell_,parameters_,title_,path_]:=Module[{s,vertexnumbers,verts,outer,inner,outerfaces,innerfaces,endfaces,innerprotoconchcenter,innerprotofaces, outerprotofaces,outerprotoconchcenter,W,T,shp,turn,Di},
{W,T,Di,turn}=parameters;

(* calculating the faces *)

vertexnumbers=Partition[Table[x,{x,Dimensions[myShell][[1]]*Dimensions[myShell][[2]]}],Dimensions[myShell][[2]]];
verts=StringRiffle[Prepend[#,"v"]]&/@Flatten[myShell,1];

outer={1,Dimensions[vertexnumbers][[2]]/2};
inner={(Dimensions[vertexnumbers][[2]]/2)+1,Dimensions[vertexnumbers][[2]]};

outerfaces=Flatten[{Table[Table[{StringRiffle[Prepend[
{vertexnumbers[[i+1,j]],vertexnumbers[[i,j+1]],vertexnumbers[[i,j]]},"f"]],StringRiffle[Prepend[{vertexnumbers[[i+1,j+1]],vertexnumbers[[i,j+1]],vertexnumbers[[i+1,j]]},"f"]]},{j,outer[[1]],outer[[2]]-1}],{i,Dimensions[vertexnumbers][[1]]-1}],Table[StringRiffle[Prepend[{vertexnumbers[[i,outer[[1]]]],vertexnumbers[[i,outer[[2]]]],vertexnumbers[[i+1,outer[[2]]]]},"f"]],{i,Dimensions[vertexnumbers][[1]]-1}],
Table[StringRiffle[Prepend[{vertexnumbers[[i,outer[[1]]]],vertexnumbers[[i+1,outer[[2]]]],vertexnumbers[[i+1,outer[[1]]]]},"f"]],{i,Dimensions[vertexnumbers][[1]]-1}]}];

innerfaces=Flatten[{Table[Table[{StringRiffle[Prepend[
{vertexnumbers[[i,j]],vertexnumbers[[i,j+1]],vertexnumbers[[i+1,j]]},"f"]],StringRiffle[Prepend[{vertexnumbers[[i+1,j]],vertexnumbers[[i,j+1]],vertexnumbers[[i+1,j+1]]},"f"]]},{j,inner[[1]],inner[[2]]-1}],{i,Dimensions[vertexnumbers][[1]]-1}],Table[StringRiffle[Prepend[{vertexnumbers[[i+1,inner[[2]]]],vertexnumbers[[i,inner[[2]]]],vertexnumbers[[i,inner[[1]]]]},"f"]],{i,Dimensions[vertexnumbers][[1]]-1}],
Table[StringRiffle[Prepend[{vertexnumbers[[i+1,inner[[1]]]],vertexnumbers[[i+1,inner[[2]]]],vertexnumbers[[i,inner[[1]]]]},"f"]],{i,Dimensions[vertexnumbers][[1]]-1}]}];

endfaces=Flatten[{Table[{StringRiffle[{"f",vertexnumbers[[-1,i+1]],vertexnumbers[[-1,i]],vertexnumbers[[-1,i+inner[[1]]-1]]}],StringRiffle[{"f",vertexnumbers[[-1,i+inner[[1]]]],vertexnumbers[[-1,i+1]],vertexnumbers[[-1,i+inner[[1]]-1]]}]},{i,outer[[2]]-1}],
StringRiffle[{"f",vertexnumbers[[-1,outer[[1]]]],vertexnumbers[[-1,outer[[2]]]],vertexnumbers[[-1,inner[[2]]]]}],
StringRiffle[{"f",vertexnumbers[[-1,outer[[1]]]],vertexnumbers[[-1,inner[[2]]]],vertexnumbers[[-1,inner[[1]]]]}]}];
innerprotoconchcenter=Mean[myShell[[1,inner[[1]];;inner[[2]]]]];
verts=Append[verts,StringRiffle[Prepend[innerprotoconchcenter,"v"]]];
innerprotofaces=Append[Table[StringRiffle[{"f",Length[verts],i+1,i}],{i,inner[[1]],inner[[2]]-1,1}],StringRiffle[{"f",Length[verts],inner[[1]],inner[[2]]}]];
outerprotoconchcenter=innerprotoconchcenter+{0,-1*Norm[innerprotoconchcenter-myShell[[1,1]]]/2,0};
verts=Append[verts,StringRiffle[Prepend[outerprotoconchcenter,"v"]]];
outerprotofaces=Append[Table[StringRiffle[{"f",Length[verts],i,i+1}],{i,outer[[1]],outer[[2]]-1,1}],StringRiffle[{"f",Length[verts],outer[[2]],outer[[1]]}]];
(* Writing the data *)
s=OpenWrite[path<>title<>".obj"];
WriteLine[s,"# Shell Mesh created by SaveShell[], P. David Polly "<>DateString[]];
WriteLine[s, "# Snails for Mathematica, Version "<>VersionNumber];
WriteLine[s,"# Parameters: W="<>ToString[W]<>" T="<>ToString[T]<>" D="<>ToString[Di]<>" turns="<>ToString[turn]];
WriteLine[s,#]&/@verts;
WriteLine[s,#]&/@outerfaces;
WriteLine[s,#]&/@innerfaces;
WriteLine[s,#]&/@endfaces;
WriteLine[s,#]&/@innerprotofaces; 
WriteLine[s,#]&/@outerprotofaces; 
WriteLine[s,"# End of File \n"];
Close[s]
];
