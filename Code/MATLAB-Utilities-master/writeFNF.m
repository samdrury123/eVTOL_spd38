function varargout = writeFNF(x,y,z,p,varargin)
%WRITEFNF
%
%   Writes a structured Cartesian grid to a PTC FNF file.
%
%   Currently only supports Quad surface cells.
%
%   Assumes data field is defined at the (i.e. cell corner structure NOT
%   cell centre).
%
%Format
%
%   "writeFNF(x,y,z,var);"
%

%default settings
nFNFrev = 3;
nElemTypes = 1;
nCSYS = 1;
nMaterials = 0;
nProperties = 0;
iiUseCSYS = 1; %sets index of CSYS to use.  N.B. CSYS #1 is defined as the absolute CSYS.
analysisName = 'Analysis1';
iiConCase = 1;

%housekeeping
currentdir = cd;

%manages variable input arguments
fname = '';
if nargin > 4
    fname = varargin{1};
end
titstr = '';
if nargin > 5
    titstr = varargin{2};
end

%establishes grid statisitcs
[nR,nC] = size(x);
nNode = nR.*nC;
nFace = (nR-1).*(nC-1);
nEdge = 2.*nFace + (nR-1) + (nC-1);
nElem = nFace;

%gets file name from user if not specified
if isempty(fname)
    [fname,pname] = uiputfile({'*.fnf','PTC FEM neutral file format (*.FNF)'},...
        'Save PTC FEM neutral file:');
    fname = fullfile(pname,fname);
    cd(currentdir) %moves back to original dir
end


%Opens file.  Return error message if there's a problem
[fid,message] = fopen(fname,'wt');
if fid < 0
    error(message);
end

%writes FNF file
writeHeader(fid);
defElements(fid);
% writeCSYS(fid); %doesn't seem to be working properly - doesn't look like
% it matters though.
writeMesh(fid,x,y,z);
writePressure(fid,p);
writeAnalysis(fid);

%ends FNF file
fprintf(fid,'%%END','');

%closes file
st = fclose(fid);
if st < 0
   warning('  Problem closing file "',fname','"'); 
end

%**************************************************************************
%NESTED-FUNCTIONS
%**************************************************************************
function writeHeader(fid)
%WRITEHEADER

%file type and date
fprintf(fid,'#PTC_FEM_NEUT %i\n',nFNFrev);
fprintf(fid,'#DATE %s\n',datestr(clock));
fprintf(fid,'\n',[]); %blank line

%file details
fprintf(fid,'%%START_SECT : %s\n','HEADER');
fprintf(fid,'%%TITLE : %s\n',titstr);
fprintf(fid,'%%STATISTICS : %i %i %i %i %i %i\n',...
    [nElemTypes,...
    nCSYS,...
    nMaterials,...
    nProperties,...
    nNode,...
    nElem]);
fprintf(fid,'%%END_SECT\n','');
fprintf(fid,'\n',[]); %blank line
    
%end nested-fucntion WRITEHEADER
end

%**************************************************************************
function defElements(fid)
%DEFELEMENTS

%starts element definition section
fprintf(fid,'%%START_SECT : %s\n','ELEM_TYPES');

%Type #1 - QUAD-SHELL definition
%------------------------------

%def
iiType.quad = 1;
fprintf(fid,'%%ELEM_TYPE %i DEF : %s %s %s %i %i %i\n',...
        iiType.quad,...
        'SHELL',... %class
        'QUAD',...  %type
        'LINEAR',...%subtype
        4, 4, 2); %no nodes, elments & faces respectively.
    
%edge connectivity
fprintf(fid,'%%ELEM_TYPE %i EDGE : %i %i %i\n',iiType.quad,1,1,2);
fprintf(fid,'%%ELEM_TYPE %i EDGE : %i %i %i\n',iiType.quad,2,2,3);
fprintf(fid,'%%ELEM_TYPE %i EDGE : %i %i %i\n',iiType.quad,3,3,4);
fprintf(fid,'%%ELEM_TYPE %i EDGE : %i %i %i\n',iiType.quad,4,4,1);

%faces
fprintf(fid,'%%ELEM_TYPE %i FACE : %i %i %i %i %i\n',iiType.quad,1,1,2,3,4);
fprintf(fid,'%%ELEM_TYPE %i FACE : %i %i %i %i %i\n',iiType.quad,1,1,4,3,2);

%ends element definition
fprintf(fid,'%%END_SECT\n','');
fprintf(fid,'\n',[]); %blank line

%end sub-fucntion DEFELEMENTS
end

%**************************************************************************
function writeCSYS(fid)
%WRITECSYS

%starts CSYS definition section
fprintf(fid,'%%START_SECT : %s\n','COORD_SYS');

%CSYS #1 - ABS coordinate system
%-------------------------------
iiCSYS = 1;
fprintf(fid,'%%COORD_SYS %i DEF : %s %s\n',iiCSYS,...
    'ABS',... %Name
    'CAR'); %Type (CARtestian, CYLindrical or SPHerical)
fprintf(fid,'%%COORD_SYS %i X_VECTOR : %f %f %f\n',iiCSYS,1,0,0);
fprintf(fid,'%%COORD_SYS %i Y_VECTOR : %f %f %f\n',iiCSYS,0,1,0);
fprintf(fid,'%%COORD_SYS %i Z_VECTOR : %f %f %f\n',iiCSYS,0,0,1);
fprintf(fid,'%%COORD_SYS %i ORIGIN : %f %f %f\n',  iiCSYS,0,0,0);

%ends CSYS definition
fprintf(fid,'%%END_SECT\n','');
fprintf(fid,'\n',[]); %blank line

%end sub-fucntion WRITECSYS
end

%**************************************************************************
function writeMesh(fid,x,y,z)
%WRITEMESH

%starts mesh definition section
fprintf(fid,'%%START_SECT : %s\n','MESH');

%writes node positions
for ii = 1:nNode
    fprintf(fid,'%%NODE %i DEF : %f %f %f %i\n',...
        ii,x(ii),y(ii),z(ii),iiUseCSYS);
end

%generates a matrix of node numbers
iiNode = [1:nNode]';
iiNode = reshape(iiNode,nR,nC);
iiCnr1 = iiNode(1:end-1, 1:end-1);
iiCnr2 = iiNode(2:end  , 1:end-1);
iiCnr3 = iiNode(2:end  , 2:end  );
iiCnr4 = iiNode(1:end-1, 2:end  );

%writes faces
for ii = 1:nFace
    fprintf(fid,'%%ELEM %i DEF : %i * * %i %i %i %i %i\n',... %N.B. material ID and property ID respectively starred out
        ii,1,iiCnr1(ii),iiCnr2(ii),iiCnr3(ii),iiCnr4(ii),iiUseCSYS);
end

%ends mesh definiton
fprintf(fid,'%%END_SECT\n','');
fprintf(fid,'\n',[]); %blank line

%end sub-fucntion WRITEMESH
end

%**************************************************************************
function writePressure(fid,p)
%WRITEPRESSURE

%default settings for a pressure field
iiLoadType = 1;
loadName = 'COEFF';
loadPlacement = 'NODE';
valueType = 'SCALAR';
iiLoad = 1;

%starts coeff field definition section
fprintf(fid,'%%START_SECT : %s\n','LOADS');

%defines load types
fprintf(fid,'%%LOAD_TYPE %i DEF : %s %s %s\n',...
    iiLoadType,loadName,loadPlacement,valueType);

%constraint case def
fprintf(fid,'%%CON_CASE %i DEF : %s\n',iiConCase,analysisName);

%defines current pressure field load
fprintf(fid,'%%LOAD %i DEF : %i %i\n',iiLoad,iiLoadType,iiConCase);

%writes pressure data
for ii = 1:nNode
   fprintf(fid,'%%LOAD %i VAL : %i %f\n',...
       iiLoad, ii, p(ii));
end

%ends coeff field definiton
fprintf(fid,'%%END_SECT\n','');
fprintf(fid,'\n',[]); %blank line

%end sub-fucntion WRITEPRESSURE
end

%**************************************************************************
function writeAnalysis(fid)
%WRITEANALYSIS

%default settings
iiSolution = 1;
solnType = 'STRUCTURAL';
solnSubType = 'STATIC';

%starts analysis definition section
fprintf(fid,'%%START_SECT : %s\n','ANALYSIS');

%analsysis definition
fprintf(fid,'%%SOLUTION %i DEF : %s %s\n',...
    iiSolution,solnType,solnSubType);
fprintf(fid,'%%SOLUTION %i CON_CASES : %i\n',iiSolution,iiConCase); 

%ends analysis definiton
fprintf(fid,'%%END_SECT\n','');
fprintf(fid,'\n',[]); %blank line

%end sub-fucntion WRITEPRESSURE
end

%**************************************************************************

%ends WRITEFNF function
end

%**************************************************************************
