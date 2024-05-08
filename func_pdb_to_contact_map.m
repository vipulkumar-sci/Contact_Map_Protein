% matlab
function native_contact_map = func_pdb_to_contact_map(pdbid, non_bonded_bondlength)
%{
    func_pdb_to_contact_map: Calculates contact map for a protein from a PDB file.
    
    Inputs:
    - pdbid: PDB file name (e.g., '2ci2.pdb').
    - non_bonded_bondlength: Threshold distance for two residues to be considered in contact.
    
    Output:
    - native_contact_map: Contact map matrix.

%}

% Opening and reading the PDB file
FileID = fopen(pdbid);
rawText = fread(FileID, inf, '*char');

% Parse lines 
splitLines = strread(rawText, '%s', 'delimiter', '\n');

% Initialize variables
numLines = length(splitLines);
PDBdata.atomType = cell(1, numLines);
PDBdata.atomNum  = cell(1, numLines);
PDBdata.atomName = cell(1, numLines);
PDBdata.resName  = cell(1, numLines);
PDBdata.resNum   = cell(1, numLines);
PDBdata.X        = cell(1, numLines);
PDBdata.Y        = cell(1, numLines);
PDBdata.Z        = cell(1, numLines);
PDBdata.comment  = cell(1, numLines);

% Loop through lines to extract relevant data
m = 1;
for n = 1:numLines
    
    thisLine = cell2mat(splitLines(n));
    
    % Check for valid line format
    if length(thisLine) > 65 && sum(isstrprop(thisLine(23:66), 'alpha')) == 0
        
        % Extract data from the line
        PDBdata.atomType{m} = thisLine(1:6);
        PDBdata.atomNum{m}  = thisLine(7:11);
        PDBdata.atomName{m} = thisLine(13:16);
        PDBdata.resName{m}  = thisLine(18:20);
        PDBdata.resNum{m}   = thisLine(23:26);
        PDBdata.X{m}        = thisLine(31:38);
        PDBdata.Y{m}        = thisLine(39:46);
        PDBdata.Z{m}        = thisLine(47:54);
        PDBdata.comment{m}  = thisLine(67:end);
        
        m = m + 1;
    end
end

% Trim excess
PDBdata.atomType(m:end) = [];
PDBdata.atomNum(m:end)  = [];
PDBdata.atomName(m:end) = [];
PDBdata.resName(m:end)  = [];
PDBdata.resNum(m:end)   = [];
PDBdata.X(m:end)        = [];
PDBdata.Y(m:end)        = [];
PDBdata.Z(m:end)        = [];
PDBdata.comment(m:end)  = [];

% Reformat data
PDBdata.atomType = strtrim(PDBdata.atomType);
PDBdata.atomNum  = str2double(PDBdata.atomNum);
PDBdata.atomName = strtrim(PDBdata.atomName);
PDBdata.resNum   = str2double(PDBdata.resNum);
PDBdata.X        = str2double(PDBdata.X);
PDBdata.Y        = str2double(PDBdata.Y);
PDBdata.Z        = str2double(PDBdata.Z);
PDBdata.comment  = strtrim(PDBdata.comment);

% Close file
fclose(FileID);

% Extract CA information for the whole chain
residue_num = [];
residue_X = [];
residue_Y = [];
residue_Z = [];
for n = 1:length(PDBdata.atomType)
    if strcmp(PDBdata.atomType{n}(1:3), 'ATO') && sum(ismember(PDBdata.atomName{n}, 'CA')) == 2
        residue_num = [residue_num; PDBdata.resNum(n)];
        residue_X = [residue_X; PDBdata.X(n)];
        residue_Y = [residue_Y; PDBdata.Y(n)];
        residue_Z = [residue_Z; PDBdata.Z(n)];
    end
end

% Initialize contact map matrix
native_contact_map = zeros(length(residue_num));

% Calculate contact map
for l = 1:length(residue_num)-1
    for m = l+2:length(residue_num)
        square_of_dist = (residue_X(l) - residue_X(m))^2 + (residue_Y(l) - residue_Y(m))^2 + (residue_Z(l) - residue_Z(m))^2;
        if square_of_dist < non_bonded_bondlength^2
            native_contact_map(l, m) = 1;
        end
    end
end

% Display contact map as heatmap
heatmap(native_contact_map);

end
