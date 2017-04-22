close all
clear all
clc

disp('Matthew C. Canver 2017')
disp('Send bugs or suggestions to Matthew_Canver AT hms DOT harvard DOT com')
disp(' ')
disp('Common PAM sequences:')
disp('S. pyogenes Cas9, PAM: NGG')
disp('S. pyogenes Cas9 variant, PAM: NGA')
disp('S. pyogenes Cas9 variant, PAM: NGCG')
disp('S. thermophilus ST1 Cas9, PAM: NNAGAA')
disp('N. meningitidis Cas9, PAM: NNNNGATT')
disp('S. aureus Cas9, PAM: NNGRRT')
disp('S. aureus Cas9 variant, PAM: NNNRRT')
disp('Acidaminococcus, Lachnospiraceae Cpf1, PAM: TTTN')
disp(' ')

%% User-specificed variables/parameters

user_set_vars_upload = textread('DNA_Striker_Variables.txt','%s');
parse_user_vars = '=';
user_set_vars = strfind(user_set_vars_upload,parse_user_vars);
index_user = find(~cellfun(@isempty,user_set_vars))+1;
list_user_vars = user_set_vars_upload(index_user);


% enter the starting and end number to correspond to the number of
% sequences to be analyzed
starting_number = str2double(list_user_vars(1,1));
end_number = str2double(list_user_vars(2,1));

% 'reference_seq_title_file' is the variable for the sequence(s) to be analyzed
% must be entered as 'filename.fasta'
reference_seq_title_file = char(list_user_vars(3,1));

% 'ref_sequence_info' is the variable for the coordinates of the sequnece(s) provided in the 'reference_seq_title_file' variable
% must be BED file format with 3 or 4 columns entered as 'filename.bed'
ref_sequence_info = char(list_user_vars(4,1));

% 'list_of_variants' is the variable containing the list of filenames
% containing the variants for each sequence(s) provided in VCF format
% entered as filename.txt. Within filename.txt, filenames should be entered
% as filenames on separate lines (do not include .VCF extension)
list_of_variants = char(list_user_vars(5,1));

PAM_sequence_list = transpose(strsplit(char(list_user_vars(6,1)),';')); % enter PAMs as {'NGG';'NGA';'NGCG'} for multiple entries or {'NGG'} for single entries
sgRNA_len_list = str2double(transpose(strsplit(char(list_user_vars(7,1)),';'))); % enter sgRNA length's as [20;23;20] for multiple entries or [20] for single entries
five_prime_PAM_list = str2double(transpose(strsplit(char(list_user_vars(8,1)),';')));; % enter 1 for 5' PAMs and 0 for 3' PAMs. enter as [0;1;0] for multiple entries or [0] for single entries
window_size = str2double(list_user_vars(9,1)); % window size for analysis
max_number_guide_ID_digits = str2double(list_user_vars(10,1)); % maximum number of digits for sgRNA ID number
coordinate_sys = char(list_user_vars(11,1)); % coordinate system chosen, enter with single quotes around it
perform_batch_output = str2double(list_user_vars(12,1)); % enter 1 to turn on batch analysis or enter 0 to turn off batch analysis
perform_variant_analysis = str2double(list_user_vars(13,1)); % enter 1 to turn on custom variant analysis or enter 0 to turn off custom variant analysis
perform_haplotype_analysis = str2double(list_user_vars(14,1)); % enter 1 to turn on haplotype analysis or enter 0 to turn off haplotype analysis
perform_plot_analysis = str2double(list_user_vars(15,1)); % enter 1 to turn on plot analysis or enter 0 to turn off plot analysis
perform_wgs_analysis = str2double(list_user_vars(16,1)); % enter 1 to turn on WGS analysis or enter 0 to turn off WGS analysis
multiple_match_analysis = str2double(list_user_vars(17,1)); % enter 1 to turn on multiple match analysis or enter 0 to turn off multiple match analysis (analyze for sgRNA matching to multiple loci within same sequence)
individual_output = str2double(list_user_vars(18,1)); % enter 1 to turn on provide individual outputs for each sequence(s) in batch analysis or enter 0 to skip individual ouputs for each sequence(s) in batch analysis
three_column_bed = str2double(list_user_vars(19,1)); % enter 1 if BED file entered as 'ref_sequence_info' has three columns: (1) chromosome, (2) start coordinate, (3) end coordinate
four_column_bed = str2double(list_user_vars(20,1)); % enter 1 if BED file entered as 'ref_sequence_info' has four columns: (1) chromosome, (2) start coordinate, (3) end coordinate, (4) other
save_figures = str2double(list_user_vars(21,1)); % enter 1 to automatically save output figures or enter 0 to skip saving figures
ks_bandwidth = str2double(list_user_vars(22,1)); % Bandwidth paramter for KS Density function (higher makes curve smoother; lower makes curve less smooth)
output_to_excel = str2double(list_user_vars(23,1)); % enter 1 to output libraries to excel or enter 0 to output analysis to a text file
haplo_guide_freq_cutoff = str2double(list_user_vars(24,1)); % enter value (0-100) to establish variant frequency cutoff for inclusion of haplotype-associated sgRNA
twenty_extra_bases = str2double(list_user_vars(25,1)); % enter 1 if included twenty extra base pairs at the ends of sequence(s) to ensure all sgRNA are obtained within interval

%% Uploading and reading reference sequence file

PAM_seq_output = strjoin(PAM_sequence_list,'_');
batch_filename = ['OUTPUT_Batch_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'];
batch_filename_gaps = ['OUTPUT_Batch_PAM=',PAM_seq_output,'_',coordinate_sys,'_gaps.xlsx'];
batch_filename_figure = ['OUTPUT_Batch_PAM=',PAM_seq_output,'_',coordinate_sys,'_Fig2.pdf'];

algorithm_type = 10000;
algorithm_type2 = 10000;
no_variant_output_num_col = 12;
haplotype_output_num_col = 8;
output_filenames_list = [];
sgRNA_gaps_full.gap_distances = [];
sgRNA_gaps_full.filenames = [];

if perform_haplotype_analysis == 1 | perform_variant_analysis == 1 | perform_wgs_analysis == 1
    vcf_filenames = textread(list_of_variants,'%s');
    if length(vcf_filenames)==1
        vcf_filenames(starting_number:end_number,1) = vcf_filenames(1,1);
    end
end

num_seq_determine = strcat(transpose(textread(reference_seq_title_file,'%c')));
num_seq_upload = length(strfind(num_seq_determine,'>'));

if num_seq_upload == 1
    [Header reference_seq1] = fastaread(reference_seq_title_file);
elseif num_seq_upload > 1
    [Header(:,1) reference_seq1(:,1)] = fastaread(reference_seq_title_file);
end

sequence_coordinates = textread(ref_sequence_info,'%s');

if four_column_bed == 1 && three_column_bed == 0
    row1 = 1;
    row2 = 2;
    row3 = 3;
    row4 = 4;
    for n = 1:(length(sequence_coordinates)/4)
        ref_seq_info1.chr(n,1) = sequence_coordinates(row1,1);
        ref_seq_info1.bpstart(n,1) = sequence_coordinates(row2,1);
        ref_seq_info1.bpend(n,1) = sequence_coordinates(row3,1);
        ref_seq_info1.identifier(n,1) = sequence_coordinates(row4,1);
        row1 = row1 + 4;
        row2 = row2 + 4;
        row3 = row3 + 4;
        row4 = row4 + 4;
    end
end

ref_seq_info.chr(starting_number:end_number,1) = ref_seq_info1.chr(starting_number:end_number);
ref_seq_info.bpstart(starting_number:end_number,1) = ref_seq_info1.bpstart(starting_number:end_number);
ref_seq_info.bpend(starting_number:end_number,1) = ref_seq_info1.bpend(starting_number:end_number);
ref_seq_info.identifier(starting_number:end_number,1) = ref_seq_info1.identifier(starting_number:end_number);

if three_column_bed == 1 && four_column_bed == 0
    row1 = 1;
    row2 = 2;
    row3 = 3;
    for n = 1:(length(sequence_coordinates)/3)
        ref_seq_info.chr(n,1) = sequence_coordinates(row1,1);
        ref_seq_info.bpstart(n,1) = sequence_coordinates(row2,1);
        ref_seq_info.bpend(n,1) = sequence_coordinates(row3,1);
        row1 = row1 + 3;
        row2 = row2 + 3;
        row3 = row3 + 3;
    end
end

if end_number > starting_number || perform_batch_output == 1
    batch_original_guides_full.guide_ID_number = [];
    batch_original_guides_full.guides_only = [];
    batch_original_guides_full.pam_only = [];
    batch_original_guides_full.pam = [];
    batch_original_guides_full.guides = [];
    batch_original_guides_full.chr_id = [];
    batch_original_guides_full.strand = [];
    batch_original_guides_full.guide_coord_start = [];
    batch_original_guides_full.guide_coord_end = [];
    batch_original_guides_full.dsb_coord = [];
    batch_original_guides_full.guide_index_start = [];
    batch_original_guides_full.guide_index_end = [];
    batch_original_guides_full.dsb_index = [];
    batch_original_guides_full.mult_match_count = [];
    batch_original_guides_full.filename = [];
end

if perform_haplotype_analysis == 1 && end_number > starting_number || perform_haplotype_analysis == 1 && perform_batch_output == 1
    batch_unique_haplotype_guides_full.guide_ID_number = [];
    batch_unique_haplotype_guides_full.guides_only = [];
    batch_unique_haplotype_guides_full.pam_only = [];
    batch_unique_haplotype_guides_full.guides = [];
    batch_unique_haplotype_guides_full.chr_id = [];
    batch_unique_haplotype_guides_full.strand = [];
    batch_unique_haplotype_guides_full.count = [];
    batch_unique_haplotype_guides_full.frequency = [];
    batch_unique_haplotype_guides_full.filename = [];
    batch_unique_haplotype_guides_full.pam_create = [];
    batch_unique_haplotype_guides_full.pam = [];
    batch_unique_haplotype_guides_full.n_mm_wt = [];
    batch_unique_haplotype_guides_full.chr_id_wt = [];
    batch_unique_haplotype_guides_full.guides_wt = [];
    batch_unique_haplotype_guides_full.guides_only_wt = [];
    batch_unique_haplotype_guides_full.pam_only_wt = [];
    batch_unique_haplotype_guides_full.strand_wt = [];
    batch_unique_haplotype_guides_full.dsb_index_wt = [];
    batch_unique_haplotype_guides_full.dsb_coord_wt = [];
    batch_unique_haplotype_guides_full.guide_index_start_wt = [];
    batch_unique_haplotype_guides_full.guide_index_end_wt = [];
    batch_unique_haplotype_guides_full.guide_coord_start_wt = [];
    batch_unique_haplotype_guides_full.guide_coord_end_wt = [];
    batch_unique_haplotype_guides_full.guide_ID_number_wt = [];
end

if multiple_match_analysis == 1 && end_number > starting_number || multiple_match_analysis == 1 && perform_batch_output == 1
    batch_multiple_match.guide_ID_number = [];
    batch_multiple_match.pam = [];
    batch_multiple_match.guides = [];
    batch_multiple_match.guides_only = [];
    batch_multiple_match.pam_only = [];
    batch_multiple_match.guides = [];
    batch_multiple_match.strand = [];
    batch_multiple_match.chr_id = [];
    batch_multiple_match.guide_coord_start = [];
    batch_multiple_match.guide_coord_end = [];
    batch_multiple_match.dsb_coord = [];
    batch_multiple_match.guide_index_start = [];
    batch_multiple_match.guide_index_end = [];
    batch_multiple_match.dsb_index = [];
    batch_multiple_match.filename = [];
end
if perform_variant_analysis == 1 && end_number > starting_number || perform_variant_analysis == 1 && perform_batch_output == 1 || perform_wgs_analysis == 1 && end_number > starting_number || perform_wgs_analysis == 1 && perform_batch_output == 1
    batch_unique_variant_guides_full.guide_ID_number = [];
    batch_unique_variant_guides_full.guides_only = [];
    batch_unique_variant_guides_full.pam_only = [];
    batch_unique_variant_guides_full.guides = [];
    batch_unique_variant_guides_full.chr_id = [];
    batch_unique_variant_guides_full.strand = [];
    batch_unique_variant_guides_full.filename = [];
    batch_unique_variant_guides_full.pam_create = [];
    batch_unique_variant_guides_full.pam = [];
    batch_unique_variant_guides_full.n_mm_wt = [];
    batch_unique_variant_guides_full.chr_id_wt = [];
    batch_unique_variant_guides_full.guides_wt = [];
    batch_unique_variant_guides_full.guides_only_wt = [];
    batch_unique_variant_guides_full.pam_only_wt = [];
    batch_unique_variant_guides_full.strand_wt = [];
    batch_unique_variant_guides_full.dsb_index_wt = [];
    batch_unique_variant_guides_full.dsb_coord_wt = [];
    batch_unique_variant_guides_full.guide_index_start_wt = [];
    batch_unique_variant_guides_full.guide_index_end_wt = [];
    batch_unique_variant_guides_full.guide_coord_start_wt = [];
    batch_unique_variant_guides_full.guide_coord_end_wt = [];
    batch_unique_variant_guides_full.guide_ID_number_wt = [];
end

for a = starting_number:end_number
    merge_original_guides_full.guide_ID_number = [];
    merge_original_guides_full.guides_only = [];
    merge_original_guides_full.pam_only = [];
    merge_original_guides_full.pam = [];
    merge_original_guides_full.guides = [];
    merge_original_guides_full.chr_id = [];
    merge_original_guides_full.strand = [];
    merge_original_guides_full.guide_coord_start = [];
    merge_original_guides_full.guide_coord_end = [];
    merge_original_guides_full.dsb_coord = [];
    merge_original_guides_full.guide_index_start = [];
    merge_original_guides_full.guide_index_end = [];
    merge_original_guides_full.dsb_index = [];
    merge_original_guides_full.mult_match_count = [];
    merge_original_guides_full.filename = [];
    if multiple_match_analysis == 1
        merge_multiple_match.guide_ID_number = [];
        merge_multiple_match.pam = [];
        merge_multiple_match.guides = [];
        merge_multiple_match.guides_only = [];
        merge_multiple_match.pam_only = [];
        merge_multiple_match.guides = [];
        merge_multiple_match.strand = [];
        merge_multiple_match.chr_id = [];
        merge_multiple_match.guide_coord_start = [];
        merge_multiple_match.guide_coord_end = [];
        merge_multiple_match.dsb_coord = [];
        merge_multiple_match.guide_index_start = [];
        merge_multiple_match.guide_index_end = [];
        merge_multiple_match.dsb_index = [];
        merge_multiple_match.filename = [];
    end
    if perform_haplotype_analysis == 1
        merge_unique_haplotype_guides_full.guide_ID_number = [];
        merge_unique_haplotype_guides_full.guides_only = [];
        merge_unique_haplotype_guides_full.pam_only = [];
        merge_unique_haplotype_guides_full.guides = [];
        merge_unique_haplotype_guides_full.chr_id = [];
        merge_unique_haplotype_guides_full.strand = [];
        merge_unique_haplotype_guides_full.count = [];
        merge_unique_haplotype_guides_full.frequency = [];
        merge_unique_haplotype_guides_full.filename = [];
        merge_unique_haplotype_guides_full.pam_create = [];
        merge_unique_haplotype_guides_full.pam = [];
        merge_unique_haplotype_guides_full.n_mm_wt = [];
        merge_unique_haplotype_guides_full.chr_id_wt = [];
        merge_unique_haplotype_guides_full.guides_wt = [];
        merge_unique_haplotype_guides_full.guides_only_wt = [];
        merge_unique_haplotype_guides_full.pam_only_wt = [];
        merge_unique_haplotype_guides_full.strand_wt = [];
        merge_unique_haplotype_guides_full.dsb_index_wt = [];
        merge_unique_haplotype_guides_full.dsb_coord_wt = [];
        merge_unique_haplotype_guides_full.guide_index_start_wt = [];
        merge_unique_haplotype_guides_full.guide_index_end_wt = [];
        merge_unique_haplotype_guides_full.guide_coord_start_wt = [];
        merge_unique_haplotype_guides_full.guide_coord_end_wt = [];
        merge_unique_haplotype_guides_full.guide_ID_number_wt = [];
    end
    
    if perform_variant_analysis == 1 || perform_wgs_analysis == 1
        merge_unique_variant_guides_full.guide_ID_number = [];
        merge_unique_variant_guides_full.guides_only = [];
        merge_unique_variant_guides_full.pam_only = [];
        merge_unique_variant_guides_full.guides = [];
        merge_unique_variant_guides_full.chr_id = [];
        merge_unique_variant_guides_full.strand = [];
        merge_unique_variant_guides_full.filename = [];
        merge_unique_variant_guides_full.pam_create = [];
        merge_unique_variant_guides_full.pam = [];
        merge_unique_variant_guides_full.n_mm_wt = [];
        merge_unique_variant_guides_full.chr_id_wt = [];
        merge_unique_variant_guides_full.guides_wt = [];
        merge_unique_variant_guides_full.guides_only_wt = [];
        merge_unique_variant_guides_full.pam_only_wt = [];
        merge_unique_variant_guides_full.strand_wt = [];
        merge_unique_variant_guides_full.dsb_index_wt = [];
        merge_unique_variant_guides_full.dsb_coord_wt = [];
        merge_unique_variant_guides_full.guide_index_start_wt = [];
        merge_unique_variant_guides_full.guide_index_end_wt = [];
        merge_unique_variant_guides_full.guide_coord_start_wt = [];
        merge_unique_variant_guides_full.guide_coord_end_wt = [];
        merge_unique_variant_guides_full.guide_ID_number_wt = [];
    end
    
    for b = 1:length(PAM_sequence_list)
        close all
        if num_seq_upload == 1;
            reference_seq_title = Header;
        elseif num_seq_upload > 1
            reference_seq_title = char(Header(a,1));
        end
        
        seq_bpstart = str2double(ref_seq_info.bpstart(a,1));
        coordinate = seq_bpstart;
        seq_bpend = str2double(ref_seq_info.bpend(a,1));
        
        disp(strcat(['Uploading and analyzing reference sequence number ',num2str(a),''],[' for'],[' ',char(PAM_sequence_list(b,1)),''],' PAM analysis...'))
        if num_seq_upload == 1;
            reference_seq = upper(reference_seq1);
        elseif num_seq_upload > 1
            reference_seq = upper(char(reference_seq1(a,1)));
        end
        
        split_reference_seq = cellstr(reference_seq')';
        
        for n = 1:length(reference_seq)
            reference_sequence.number(n,1) = n;
            reference_sequence.coord(n,1) = coordinate;
            coordinate = coordinate+1;
            reference_sequence.nucleotide(n,1) = split_reference_seq(n);
        end
        
        %% Set PAM Sequence and Identify All Possible sgRNA
        
        PAM_sequence1 = char(PAM_sequence_list(b,1));
        sgRNA_len = sgRNA_len_list(b,1);
        five_prime_PAM = five_prime_PAM_list(b,1);
        
        PAM_len = length(PAM_sequence1);
        N_seq_in_PAM = length(strfind(PAM_sequence1,'N'));
        if five_prime_PAM == 0
            index_PAM = N_seq_in_PAM+1;
        end
        if five_prime_PAM == 1
            index_PAM = 1;
        end
        
        if five_prime_PAM == 0
            if length(cellstr(PAM_sequence1)) == 1
                PAM_sequence = PAM_sequence1(index_PAM:PAM_len);
            elseif length(cellstr(PAM_sequence1)) > 1
                PAM_len = length(PAM_sequence1{1,1});
                N_seq_in_PAM = length(strfind(PAM_sequence1{1,1},'N'));
                index_PAM = N_seq_in_PAM+1;
                for n = 1:length(cellstr(PAM_sequence1))
                    PAM_seq_temp = char(PAM_sequence1(n,1));
                    PAM_sequence{n,1} = PAM_seq_temp(index_PAM:PAM_len);
                end
            end
        end
        
        if five_prime_PAM == 1
            if length(cellstr(PAM_sequence1)) == 1
                PAM_sequence = PAM_sequence1(index_PAM:(PAM_len-N_seq_in_PAM));
            elseif length(cellstr(PAM_sequence1)) > 1
                PAM_len = length(PAM_sequence1{1,1});
                N_seq_in_PAM = length(strfind(PAM_sequence1{1,1},'N'));
                index_PAM = 1;
                for n = 1:length(cellstr(PAM_sequence1))
                    PAM_seq_temp = char(PAM_sequence1(n,1));
                    PAM_sequence{n,1} = PAM_seq_temp(index_PAM:(PAM_len-N_seq_in_PAM));
                end
            end
        end
        
        if sum(ismember(PAM_sequence,'R')) > 0
            R_seq_in_PAM(:,1) = strfind(PAM_sequence,'R');
        end
        
        if exist('R_seq_in_PAM','var')
            row = 1;
            PAM_sequence_with_R = PAM_sequence;
            PAM_sequence_with_R(R_seq_in_PAM) = 'A';
            PAM_sequence_1{row,1} = PAM_sequence_with_R;
            PAM_sequence_1{row,1} = PAM_sequence_with_R;
            row = row + 1;
            PAM_sequence_with_R(R_seq_in_PAM) = 'G';
            PAM_sequence_1{row,1} = PAM_sequence_with_R;
            row = row + 1;
            if length(R_seq_in_PAM)>1
                for n = 1:length(R_seq_in_PAM)
                    PAM_sequence_with_R = PAM_sequence_1{2,1};
                    PAM_sequence_with_R(R_seq_in_PAM(n,1)) = 'A';
                    PAM_sequence_1{row,1} = PAM_sequence_with_R;
                    row = row + 1;
                end
                PAM_sequence = PAM_sequence_1;
            else PAM_sequence = PAM_sequence_1;
            end
        end
        
        reference_sequence_rc_full = seqrcomplement(strjoin(transpose(reference_sequence.nucleotide),''));
        split_reference_seq_rc = cellstr(reference_sequence_rc_full')';
        reference_sequence_rc.coord = reference_sequence.coord;
        reference_sequence_rc.number_for = reference_sequence.number;
        reference_sequence_rc.number_rev = sort(reference_sequence.number,'descend');
        reference_sequence_rc.nucleotide = transpose(split_reference_seq_rc);
        
        guide_index1 = [];
        
        if sum(ismember(PAM_sequence1,'R')) == 0 & length(cellstr(PAM_sequence)) == 1
            PAM_sequence_cell = {PAM_sequence};
        else PAM_sequence_cell = PAM_sequence;
        end
        
        %% Custom variant analysis
        
        if perform_variant_analysis == 1 || perform_wgs_analysis == 1
            if perform_variant_analysis == 1 && perform_wgs_analysis ==  0;
                num_vcf_columns = 8;
                for n = 1:num_vcf_columns
                    eval(['row',num2str(n),'=',num2str(n),';']);
                    if n == num_vcf_columns
                        eval(['row',num2str(n+1),'=',num2str(n+1),';'])
                    end
                end
                temp_vcf_name  = [vcf_filenames(a,1),'.vcf'];
                vcf_file4 = textread(char(strjoin(temp_vcf_name,'')),'%s');
                vcf_file3 = strfind(vcf_file4,'#CHROM');
                vcf_file2 = find(~cellfun(@isempty,vcf_file3));
                vcf_file1 = vcf_file4(vcf_file2:length(vcf_file4));
                vcf_file = vcf_file1(num_vcf_columns+1:length(vcf_file1));
                
                for n = 1:length(vcf_file)
                    if row1 <= length(vcf_file)
                        eval(['vcf_snp.chr_id(n,1) = str2double(vcf_file(row1,1));']);
                        eval(['vcf_snp.bpstart(n,1) = str2double(vcf_file(row2,1));']);
                        eval(['vcf_snp.rs_id(n,1) = vcf_file(row3,1);']);
                        eval(['vcf_snp.ref_snp(n,1) = vcf_file(row4,1);']);
                        eval(['vcf_snp.alt_snp(n,1) = vcf_file(row5,1);']);
                        eval(['vcf_snp.qual(n,1) = vcf_file(row6,1);']);
                        eval(['vcf_snp.filter(n,1) = vcf_file(row7,1);']);
                        if strfind(char(vcf_file(row8,1)),'GMAF=')>0
                            gmaf_temp1 = strsplit(char(vcf_file(row8,1)),'GMAF=');
                            gmaf_temp2 = strsplit(char(gmaf_temp1(1,2)),';VC');
                            eval(['vcf_snp.gmaf(n,1) = str2double(gmaf_temp2(1,1));']);
                        elseif isempty(strfind(char(vcf_file(row8,1)),'GMAF=')) == 1;
                            vcf_snp.gmaf(n,1) = 0;
                        end
                        eval(['vcf_snp.info(n,1) = vcf_file(row9,1);']);
                        row1 = row1+num_vcf_columns;
                        row2 = row2+num_vcf_columns;
                        row3 = row3+num_vcf_columns;
                        row4 = row4+num_vcf_columns;
                        row5 = row5+num_vcf_columns;
                        row6 = row6+num_vcf_columns;
                        row7 = row7+num_vcf_columns;
                        row8 = row8+num_vcf_columns;
                    elseif row1>length(vcf_file)
                        break
                    end
                end
            end
            
            %% WGS Seq
            
            if perform_wgs_analysis == 1 && perform_variant_analysis == 0
                num_vcf_columns = 10;
                for n = 1:num_vcf_columns
                    eval(['row',num2str(n),'=',num2str(n),';']);
                    if n == num_vcf_columns
                        eval(['row',num2str(n+1),'=',num2str(n+1),';'])
                    end
                end
                temp_vcf_name  = [vcf_filenames(a,1),'.vcf'];
                vcf_file4 = textread(char(strjoin(temp_vcf_name,'')),'%s');
                vcf_file3 = strfind(vcf_file4,'#CHROM');
                vcf_file2 = find(~cellfun(@isempty,vcf_file3));
                vcf_file1 = vcf_file4(vcf_file2:length(vcf_file4));
                vcf_file = vcf_file1(num_vcf_columns+1:length(vcf_file1));
                
                for n = 1:length(vcf_file)
                    if row1 <= length(vcf_file)
                        eval(['vcf_snp.chr_id(n,1) = str2double(vcf_file(row1,1));']);
                        eval(['vcf_snp.bpstart(n,1) = str2double(vcf_file(row2,1));']);
                        eval(['vcf_snp.rs_id(n,1) = vcf_file(row3,1);']);
                        eval(['vcf_snp.ref_snp(n,1) = vcf_file(row4,1);']);
                        eval(['vcf_snp.alt_snp(n,1) = vcf_file(row5,1);']);
                        eval(['vcf_snp.qual(n,1) = vcf_file(row6,1);']);
                        eval(['vcf_snp.filter(n,1) = vcf_file(row7,1);']);
                        eval(['vcf_snp.info(n,1) = vcf_file(row8,1);']);
                        eval(['vcf_snp.format(n,1) = vcf_file(row9,1);']);
                        eval(['vcf_snp.wgs_id(n,1) = vcf_file(row10,1);']);
                        row1 = row1+num_vcf_columns;
                        row2 = row2+num_vcf_columns;
                        row3 = row3+num_vcf_columns;
                        row4 = row4+num_vcf_columns;
                        row5 = row5+num_vcf_columns;
                        row6 = row6+num_vcf_columns;
                        row7 = row7+num_vcf_columns;
                        row8 = row8+num_vcf_columns;
                        row9 = row9+num_vcf_columns;
                        row10 = row10+num_vcf_columns;
                    elseif row1>length(vcf_file)
                        break
                    end
                end
            end
            
            %% Sorting variants from lowest to highest coordinate
            
            [sorted_bpstart index_sorting] = sort(vcf_snp.bpstart,1);
            sorted_vcf_snp.chr_id = vcf_snp.chr_id(index_sorting);
            sorted_vcf_snp.bpstart = vcf_snp.bpstart(index_sorting);
            sorted_vcf_snp.rs_id = vcf_snp.rs_id(index_sorting);
            sorted_vcf_snp.ref_snp = vcf_snp.ref_snp(index_sorting);
            sorted_vcf_snp.alt_snp = vcf_snp.alt_snp(index_sorting);
            sorted_vcf_snp.qual = vcf_snp.qual(index_sorting);
            sorted_vcf_snp.filter = vcf_snp.filter(index_sorting);
            sorted_vcf_snp.info = vcf_snp.info(index_sorting);
            
            wgs_custom_var = 1;
            relevant_snps_index = find(sorted_vcf_snp.bpstart>=min(reference_sequence.coord) & sorted_vcf_snp.bpstart<=max(reference_sequence.coord));
            sorted_vcf_snp.chr_id = sorted_vcf_snp.chr_id(relevant_snps_index);
            sorted_vcf_snp.bpstart = sorted_vcf_snp.bpstart(relevant_snps_index);
            sorted_vcf_snp.rs_id = sorted_vcf_snp.rs_id(relevant_snps_index);
            sorted_vcf_snp.ref_snp = sorted_vcf_snp.ref_snp(relevant_snps_index);
            sorted_vcf_snp.alt_snp = sorted_vcf_snp.alt_snp(relevant_snps_index);
            sorted_vcf_snp.qual = sorted_vcf_snp.qual(relevant_snps_index);
            sorted_vcf_snp.filter = sorted_vcf_snp.filter(relevant_snps_index);
            sorted_vcf_snp.info = sorted_vcf_snp.info(relevant_snps_index);
            if length(relevant_snps_index)==0
                wgs_custom_var = 0;
            end
            
            comma_alt = strfind(sorted_vcf_snp.alt_snp,',');
            comma_alt_index = find(~cellfun(@isempty,comma_alt));
            
            if length(comma_alt_index)>0
            alt_row = 1;
            for r = 1:length(comma_alt_index)
                split_var_temp = transpose(strsplit(char(sorted_vcf_snp.alt_snp(comma_alt_index(r,1),1)),','));
                sorted_vcf_snp.alt_snp(comma_alt_index(r,1),1) = split_var_temp(1,1);
                alt_sorted_vcf_snp.alt_snp(alt_row,1) = split_var_temp(2,1);
                alt_sorted_vcf_snp.chr_id(alt_row,1) = sorted_vcf_snp.chr_id(comma_alt_index(r,1));
                alt_sorted_vcf_snp.bpstart(alt_row,1) = sorted_vcf_snp.bpstart(comma_alt_index(r,1));
                alt_sorted_vcf_snp.rs_id(alt_row,1) = sorted_vcf_snp.rs_id(comma_alt_index(r,1));
                alt_sorted_vcf_snp.ref_snp(alt_row,1) = sorted_vcf_snp.ref_snp(comma_alt_index(r,1));
                alt_sorted_vcf_snp.qual(alt_row,1) = sorted_vcf_snp.qual(comma_alt_index(r,1));
                alt_sorted_vcf_snp.filter(alt_row,1) = sorted_vcf_snp.filter(comma_alt_index(r,1));
                alt_sorted_vcf_snp.info(alt_row,1) = sorted_vcf_snp.info(comma_alt_index(r,1));
                alt_row = alt_row + 1;
            end
            
            sorted_vcf_snp.chr_id = [sorted_vcf_snp.chr_id;alt_sorted_vcf_snp.chr_id];
            sorted_vcf_snp.bpstart = [sorted_vcf_snp.bpstart;alt_sorted_vcf_snp.bpstart];
            sorted_vcf_snp.rs_id =[sorted_vcf_snp.rs_id;alt_sorted_vcf_snp.rs_id];
            sorted_vcf_snp.ref_snp = [sorted_vcf_snp.ref_snp;alt_sorted_vcf_snp.ref_snp];
            sorted_vcf_snp.alt_snp = [sorted_vcf_snp.alt_snp;alt_sorted_vcf_snp.alt_snp];
            sorted_vcf_snp.qual = [sorted_vcf_snp.qual;alt_sorted_vcf_snp.qual];
            sorted_vcf_snp.filter = [sorted_vcf_snp.filter;alt_sorted_vcf_snp.filter];
            sorted_vcf_snp.info = [sorted_vcf_snp.info;alt_sorted_vcf_snp.info];
            end
            
            %% Classifying variants as substitutions, deletions, or insertions
            nucleotide_list = ['AGCT'];
            row_sub = 1;
            row_del = 1;
            row_ins = 1;
            if wgs_custom_var == 1;
                for n = 1:length(sorted_vcf_snp.bpstart)
                    temp1 = strfind(nucleotide_list,char(sorted_vcf_snp.ref_snp(n,1)));
                    temp2 = strfind(nucleotide_list,char(sorted_vcf_snp.alt_snp(n,1)));
                    if length(find(ismember(nucleotide_list,char(sorted_vcf_snp.ref_snp(n,1)))>0))>0 & length(find(ismember(nucleotide_list,char(sorted_vcf_snp.alt_snp(n,1)))>0))>0
                        substitutions.chr_id(row_sub,1) = sorted_vcf_snp.chr_id(n,1);
                        substitutions.bpstart(row_sub,1) = sorted_vcf_snp.bpstart(n,1);
                        substitutions.bpend(row_sub,1) = sorted_vcf_snp.bpstart(n,1);
                        substitutions.rs_id(row_sub,1) = sorted_vcf_snp.rs_id(n,1);
                        substitutions.ref_snp(row_sub,1) = sorted_vcf_snp.ref_snp(n,1);
                        substitutions.alt_snp(row_sub,1) = sorted_vcf_snp.alt_snp(n,1);
                        substitutions.qual(row_sub,1) = sorted_vcf_snp.qual(n,1);
                        substitutions.filter(row_sub,1) = sorted_vcf_snp.filter(n,1);
                        substitutions.info(row_sub,1) = sorted_vcf_snp.info(n,1);
                        substitutions.type(row_sub,1) = {'Substitution'};
                        substitutions.identifier(row_sub,1) = 1;
                        row_sub=row_sub+1;
                    elseif length(find(ismember('-',char(sorted_vcf_snp.alt_snp(n,1)))>0))>0
                        deletions.chr_id(row_del,1) = sorted_vcf_snp.chr_id(n,1);
                        deletions.bpstart(row_del,1) = sorted_vcf_snp.bpstart(n,1);
                        deletions.rs_id(row_del,1) = sorted_vcf_snp.rs_id(n,1);
                        deletions.ref_snp(row_del,1) = sorted_vcf_snp.ref_snp(n,1);
                        deletions.alt_snp(row_del,1) = sorted_vcf_snp.alt_snp(n,1);
                        deletions.qual(row_del,1) = sorted_vcf_snp.qual(n,1);
                        deletions.filter(row_del,1) = sorted_vcf_snp.filter(n,1);
                        deletions.info(row_del,1) = sorted_vcf_snp.info(n,1);
                        deletions.type(row_del,1) = {'Deletion'};
                        deletions.identifier(row_del,1) = 2;
                        deletions.del_len(row_del,1) = size(char(sorted_vcf_snp.ref_snp(n,1)),2);
                        deletions.bpend(row_del,1) = deletions.bpstart(row_del,1)+deletions.del_len(row_del,1)-1;
                        row_del=row_del+1;
                    elseif length(find(ismember('-',char(sorted_vcf_snp.ref_snp(n,1)))>0))>0
                        insertions.chr_id(row_ins,1) = sorted_vcf_snp.chr_id(n,1);
                        insertions.bpstart(row_ins,1) = sorted_vcf_snp.bpstart(n,1);
                        insertions.rs_id(row_ins,1) = sorted_vcf_snp.rs_id(n,1);
                        insertions.ref_snp(row_ins,1) = sorted_vcf_snp.ref_snp(n,1);
                        insertions.alt_snp(row_ins,1) = sorted_vcf_snp.alt_snp(n,1);
                        insertions.qual(row_ins,1) = sorted_vcf_snp.qual(n,1);
                        insertions.filter(row_ins,1) = sorted_vcf_snp.filter(n,1);
                        insertions.info(row_ins,1) = sorted_vcf_snp.info(n,1);
                        insertions.type(row_ins,1) = {'Insertion'};
                        insertions.identifier(row_ins,1) = 3;
                        insertions.ins_len(row_ins,1) = size(char(sorted_vcf_snp.alt_snp(n,1)),2);
                        insertions.bpend(row_ins,1) = insertions.bpstart(row_ins,1)+insertions.ins_len(row_ins,1)-1;
                        row_ins=row_ins+1;
                    end
                end
            end
            %Substitutions
            if exist('substitutions','var') == 1;
                for n = 1:length(substitutions.bpstart)
                    substitutions.seq_index(n,1) = find(substitutions.bpstart(n,1)==reference_sequence.coord);
                end
            end
            
            %Deletions
            if exist('deletions','var') == 1;
                for n = 1:length(deletions.bpstart)
                    deletions.seq_index_start(n,1) = find(deletions.bpstart(n,1) == reference_sequence.coord);
                    deletions.seq_index_end(n,1) = find(deletions.bpend(n,1) == reference_sequence.coord);
                end
            end
            
            %Insertions
            if exist('insertions','var') == 1;
                for n = 1:length(insertions.bpstart)
                    insertions.seq_index_start(n,1) = find(insertions.bpstart(n,1) == reference_sequence.coord);
                    insertions.seq_index_end(n,1) = find(insertions.bpend(n,1) == reference_sequence.coord);
                end
            end
            
            test_variants.bpstart = [];
            test_variants.bpend = [];
            test_variants.rs_id = [];
            test_variants.ref_snp = [];
            test_variants.alt_snp = [];
            test_variants.type = [];
            test_variants.identifier = [];
            test_variants.variant_index_start = [];
            test_variants.variant_index_end = [];
            test_variants.gmaf = [];
            
            if exist('substitutions','var') == 1
                test_variants.bpstart = [test_variants.bpstart;substitutions.bpstart];
                test_variants.bpend = [test_variants.bpend;substitutions.bpend];
                test_variants.rs_id = [test_variants.rs_id;substitutions.rs_id];
                test_variants.type = [test_variants.type;substitutions.type];
                test_variants.ref_snp = [test_variants.ref_snp;substitutions.ref_snp];
                test_variants.alt_snp = [test_variants.alt_snp;substitutions.alt_snp];
                test_variants.identifier = [test_variants.identifier;substitutions.identifier];
                test_variants.variant_index_start = [test_variants.variant_index_start;substitutions.seq_index];
                test_variants.variant_index_end = [test_variants.variant_index_end;substitutions.seq_index];
            end
            if exist('deletions','var') == 1
                test_variants.bpstart = [test_variants.bpstart;deletions.bpstart];
                test_variants.bpend = [test_variants.bpend;deletions.bpend];
                test_variants.rs_id = [test_variants.rs_id;deletions.rs_id];
                test_variants.type = [test_variants.type;deletions.type];
                test_variants.ref_snp = [test_variants.ref_snp;deletions.ref_snp];
                test_variants.alt_snp = [test_variants.alt_snp;deletions.alt_snp];
                test_variants.identifier = [test_variants.identifier;deletions.identifier];
                test_variants.variant_index_start = [test_variants.variant_index_start;deletions.seq_index_start];
                test_variants.variant_index_end = [test_variants.variant_index_end;deletions.seq_index_end];
            end
            if exist('insertions','var') == 1
                test_variants.bpstart = [test_variants.bpstart;insertions.bpstart];
                test_variants.bpend = [test_variants.bpend;insertions.bpend];
                test_variants.rs_id = [test_variants.rs_id;insertions.rs_id];
                test_variants.type = [test_variants.type;insertions.type];
                test_variants.ref_snp = [test_variants.ref_snp;insertions.ref_snp];
                test_variants.alt_snp = [test_variants.alt_snp;insertions.alt_snp];
                test_variants.identifier = [test_variants.identifier;insertions.identifier];
                test_variants.variant_index_start = [test_variants.variant_index_start;insertions.seq_index_start];
                test_variants.variant_index_end = [test_variants.variant_index_end;insertions.seq_index_end];
            end
            
            [test_variants.bpstart index_test_variants_sorting] = sort(test_variants.bpstart,1);
            test_variants.bpend = test_variants.bpend(index_test_variants_sorting);
            test_variants.rs_id = test_variants.rs_id(index_test_variants_sorting);
            test_variants.ref_snp = test_variants.ref_snp(index_test_variants_sorting);
            test_variants.alt_snp = test_variants.alt_snp(index_test_variants_sorting);
            test_variants.identifier = test_variants.identifier(index_test_variants_sorting);
            test_variants.variant_index_start = test_variants.variant_index_start(index_test_variants_sorting);
            test_variants.variant_index_end = test_variants.variant_index_end(index_test_variants_sorting);
            test_variants.variant_len = test_variants.bpend-test_variants.bpstart;
            
            test_variants.chr_id(1:length(test_variants.bpstart),1) = ref_seq_info.chr(a,1);
            
            for n = 1:length(test_variants.variant_len)
                if test_variants.identifier(n,1) == 2
                    test_variants.variant_len_full(n,1) = -1*(test_variants.variant_len(n,1)+1);
                elseif test_variants.identifier(n,1) == 3
                    test_variants.variant_len_full(n,1) = test_variants.variant_len(n,1)+1;
                elseif test_variants.identifier(n,1) == 1
                    test_variants.variant_len_full(n,1) = 0;
                end
            end
            
        end
        
        %% Sliding Window Approach
        window_start = 1;
        window_end = window_size+1;
        num_windows = 2*floor(length(reference_seq)/window_size)+1;
        guide_num = seq_bpstart;
        guide_num_index = 1;
        semicolon = ';';
        
        for n = 1:num_windows
            clearvars variants_in_window window_sequence window_sequence_full window_sequence_sub window_sequence_del window_sequence_ins near_variants
            window_delta = 0.5*(n-1)*window_size;
            window_start = 1+window_delta;
            window_end = (window_size+1)+window_delta;
            if window_end > length(reference_seq)
                window_end = length(reference_seq);
            end
            window_sequence_full = strjoin(transpose(reference_sequence.nucleotide(window_start:window_end)),'');
            window_sequence.nucleotide = transpose(cellstr(window_sequence_full')');
            window_sequence.index = reference_sequence.number(window_start:window_end);
            window_sequence.coord = reference_sequence.coord(window_start:window_end);
            guide_num = seq_bpstart+window_delta;
            grade_num_index = guide_num_index+window_delta;
            if perform_variant_analysis == 1 && length(test_variants.variant_index_start)>0 || perform_wgs_analysis == 1 && length(test_variants.variant_index_start)>0
                row = 1;
                for m = 1:length(test_variants.variant_index_start)
                    idx1(m,1) = (test_variants.variant_index_start(m)-window_start).*(window_end-test_variants.variant_index_end(m));
                end
                idx = find(idx1>0);
                if length(idx)>0
                    variants_in_window.ref_snp = test_variants.ref_snp(idx,1);
                    variants_in_window.alt_snp = test_variants.alt_snp(idx,1);
                    variants_in_window.variant_index_start = test_variants.variant_index_start(idx,1);
                    variants_in_window.variant_index_end = test_variants.variant_index_end(idx,1);
                    variants_in_window.variant_coord_start = test_variants.bpstart(idx,1);
                    variants_in_window.variant_coord_end = test_variants.bpend(idx,1);
                    variants_in_window.rs_id = test_variants.rs_id(idx,1);
                    variants_in_window.identifier = test_variants.identifier(idx,1);
                    variants_in_window.variant_len = test_variants.variant_len(idx,1);
                end
                
                if exist('variants_in_window','var') == 1;
                    [variants_in_window.variant_coord_start sort_index] = sort(variants_in_window.variant_coord_start,1);
                    variants_in_window.ref_snp = variants_in_window.ref_snp(sort_index);
                    variants_in_window.alt_snp = variants_in_window.alt_snp(sort_index);
                    variants_in_window.variant_index_start = variants_in_window.variant_index_start(sort_index);
                    variants_in_window.variant_index_end = variants_in_window.variant_index_end(sort_index);
                    variants_in_window.variant_coord_start = variants_in_window.variant_coord_start(sort_index);
                    variants_in_window.variant_coord_end = variants_in_window.variant_coord_end(sort_index);
                    variants_in_window.rs_id = variants_in_window.rs_id(sort_index);
                    variants_in_window.identifier = variants_in_window.identifier(sort_index);
                    variants_in_window.variant_len = variants_in_window.variant_len(sort_index);
                end
            end
            if exist('variants_in_window','var') == 1;
                num_variants_window = length(variants_in_window.ref_snp);
            else num_variants_window = 0;
            end
            num_variants_per_window(n,1) = num_variants_window;
            row_w = 1;
            eval(['window_sequence_var',num2str(n),'.window{row_w,1} = window_sequence_full;']);
            eval(['window_sequence_var',num2str(n),'.window_start(row_w,1) = window_start;']);
            eval(['window_sequence_var',num2str(n),'.window_end(row_w,1) = window_end;']);
            row_w = row_w + 1;
            if perform_variant_analysis == 1 && perform_wgs_analysis == 0
                if num_variants_window > 0;
                    for j = 1:length(variants_in_window.ref_snp);
                        if variants_in_window.identifier(j,1) == 1 & strcmp(window_sequence.nucleotide(find(variants_in_window.variant_index_start(j,1)==window_sequence.index)),variants_in_window.ref_snp(j,1)) == 1;
                            window_sequence_sub = window_sequence.nucleotide(:,1);
                            window_sequence_sub(find(variants_in_window.variant_index_start(j,1)==window_sequence.index),1) = variants_in_window.alt_snp(j,1);
                            eval(['window_sequence_var',num2str(n),'.window{row_w,1} = strjoin(transpose(window_sequence_sub),'''');']);
                            eval(['window_sequence_var',num2str(n),'.ref_snp(row_w,1) = variants_in_window.ref_snp(j,1);']);
                            eval(['window_sequence_var',num2str(n),'.alt_snp(row_w,1) = variants_in_window.alt_snp(j,1);']);
                            eval(['window_sequence_var',num2str(n),'.from_seq(row_w,1) = window_sequence.nucleotide(find(variants_in_window.variant_index_start(j,1)==window_sequence.index));']);
                            eval(['window_sequence_var',num2str(n),'.variant_index_start(row_w,1) = variants_in_window.variant_index_start(j,1);']);
                            eval(['window_sequence_var',num2str(n),'.variant_index_end(row_w,1) = variants_in_window.variant_index_end(j,1);']);
                            eval(['window_sequence_var',num2str(n),'.variant_coord_start(row_w,1) = variants_in_window.variant_coord_start(j,1);']);
                            eval(['window_sequence_var',num2str(n),'.variant_coord_end(row_w,1) = variants_in_window.variant_coord_end(j,1);']);
                            eval(['window_sequence_var',num2str(n),'.rs_id(row_w,1) = variants_in_window.rs_id(j,1);']);
                            eval(['window_sequence_var',num2str(n),'.identifier(row_w,1) = variants_in_window.identifier(j,1);']);
                            eval(['window_sequence_var',num2str(n),'.variant_len(row_w,1) = variants_in_window.variant_len(j,1);']);
                            eval(['window_sequence_var',num2str(n),'.type{row_w,1} = ''Substitution'';']);
                            eval(['window_sequence_var',num2str(n),'.num_var{row_w,1} = ''Single Variant'';']);
                            row_w = row_w + 1;
                        elseif variants_in_window.identifier(j,1) == 2 & strcmp(strjoin(transpose(window_sequence.nucleotide(find(variants_in_window.variant_index_start(j,1)==window_sequence.index):find(variants_in_window.variant_index_end(j,1)==window_sequence.index))),''),variants_in_window.ref_snp(j,1)) == 1;
                            eval(['window_sequence_var',num2str(n),'.ref_snp(row_w,1) = variants_in_window.ref_snp(j,1);']);
                            eval(['window_sequence_var',num2str(n),'.alt_snp(row_w,1) = variants_in_window.alt_snp(j,1);']);
                            eval(['window_sequence_var',num2str(n),'.from_seq(row_w,1) = window_sequence.nucleotide(find(variants_in_window.variant_index_start(j,1)==window_sequence.index));']);
                            eval(['window_sequence_var',num2str(n),'.variant_index_start(row_w,1) = variants_in_window.variant_index_start(j,1);']);
                            eval(['window_sequence_var',num2str(n),'.variant_index_end(row_w,1) = variants_in_window.variant_index_end(j,1);']);
                            eval(['window_sequence_var',num2str(n),'.variant_coord_start(row_w,1) = variants_in_window.variant_coord_start(j,1);']);
                            eval(['window_sequence_var',num2str(n),'.variant_coord_end(row_w,1) = variants_in_window.variant_coord_end(j,1);']);
                            eval(['window_sequence_var',num2str(n),'.rs_id(row_w,1) = variants_in_window.rs_id(j,1);']);
                            eval(['window_sequence_var',num2str(n),'.identifier(row_w,1) = variants_in_window.identifier(j,1);']);
                            eval(['window_sequence_var',num2str(n),'.variant_len(row_w,1) = variants_in_window.variant_len(j,1);']);
                            eval(['window_sequence_var',num2str(n),'.type{row_w,1} = ''Deletion'';']);
                            eval(['window_sequence_var',num2str(n),'.num_var{row_w,1} = ''Single Variant'';']);
                            window_sequence_del = window_sequence.nucleotide(:,1);
                            window_sequence_del(find(variants_in_window.variant_index_start(j,1)==window_sequence.index):find(variants_in_window.variant_index_end(j,1)==window_sequence.index)) = {'D'};
                            ident_del_sites = strfind(window_sequence_del,'D');
                            ident_del_sites = find(~cellfun(@isempty,ident_del_sites));
                            window_sequence_del(ident_del_sites) = [];
                            eval(['window_sequence_var',num2str(n),'.window{row_w,1} = strjoin(transpose(window_sequence_del),'''');']);
                            row_w = row_w + 1;
                        elseif variants_in_window.identifier(j,1) == 3
                            eval(['window_sequence_var',num2str(n),'.ref_snp(row_w,1) = variants_in_window.ref_snp(j,1);']);
                            eval(['window_sequence_var',num2str(n),'.alt_snp(row_w,1) = variants_in_window.alt_snp(j,1);']);
                            eval(['window_sequence_var',num2str(n),'.from_seq(row_w,1) = window_sequence.nucleotide(find(variants_in_window.variant_index_start(j,1)==window_sequence.index));']);
                            eval(['window_sequence_var',num2str(n),'.variant_index_start(row_w,1) = variants_in_window.variant_index_start(j,1);']);
                            eval(['window_sequence_var',num2str(n),'.variant_index_end(row_w,1) = variants_in_window.variant_index_end(j,1);']);
                            eval(['window_sequence_var',num2str(n),'.variant_coord_start(row_w,1) = variants_in_window.variant_coord_start(j,1);']);
                            eval(['window_sequence_var',num2str(n),'.variant_coord_end(row_w,1) = variants_in_window.variant_coord_end(j,1);']);
                            eval(['window_sequence_var',num2str(n),'.rs_id(row_w,1) = variants_in_window.rs_id(j,1);']);
                            eval(['window_sequence_var',num2str(n),'.identifier(row_w,1) = variants_in_window.identifier(j,1);']);
                            eval(['window_sequence_var',num2str(n),'.variant_len(row_w,1) = variants_in_window.variant_len(j,1);']);
                            eval(['window_sequence_var',num2str(n),'.type{row_w,1} = ''Insertion'';']);
                            eval(['window_sequence_var',num2str(n),'.num_var{row_w,1} = ''Single Variant'';']);
                            window_sequence_ins = window_sequence.nucleotide(:,1);
                            ins_length_increase = 0;
                            split_ins = transpose(cellstr(char(variants_in_window.alt_snp(j,1))')');
                            ins_length = length(split_ins);
                            index_ins_original = find(variants_in_window.variant_index_start(j,1)==window_sequence.index);
                            index_ins = find(variants_in_window.variant_index_start(j,1)==window_sequence.index)+ins_length_increase;
                            window_sequence_ins_final = window_sequence_ins(1:(index_ins-1));
                            window_sequence_ins_final((index_ins+ins_length):(length(window_sequence_ins)+ins_length)) = window_sequence_ins(index_ins:length(window_sequence_ins));
                            ins_length_increase = ins_length_increase + 1;
                            index_replace = index_ins;
                            ins_count = 0;
                            for i = 1:ins_length
                                window_sequence_ins_final(index_replace+ins_count) = split_ins(i,1);
                                ins_count = ins_count + 1;
                            end
                            eval(['window_sequence_var',num2str(n),'.window{row_w,1} = strjoin(transpose(window_sequence_ins_final),'''');']);
                            row_w = row_w + 1;
                        end
                    end
                end
            end
            if perform_wgs_analysis == 1 && perform_variant_analysis == 0
                if num_variants_window > 0;
                    window_sequence_mod = window_sequence.nucleotide(:,1);
                    for j = 1:length(variants_in_window.ref_snp);
                        if variants_in_window.identifier(j,1) == 1 & strcmp(window_sequence.nucleotide(find(variants_in_window.variant_index_start(j,1)==window_sequence.index)),variants_in_window.ref_snp(j,1)) == 1;
                            window_sequence_mod(find(variants_in_window.variant_index_start(j,1)==window_sequence.index),1) = variants_in_window.alt_snp(j,1);
                            eval(['window_sequence_var',num2str(n),'.window{row_w,1} = strjoin(transpose(window_sequence_mod),'''');']);
                            eval(['window_sequence_var',num2str(n),'.ref_snp(row_w,1) = variants_in_window.ref_snp(j,1);']);
                            eval(['window_sequence_var',num2str(n),'.alt_snp(row_w,1) = variants_in_window.alt_snp(j,1);']);
                            eval(['window_sequence_var',num2str(n),'.from_seq(row_w,1) = window_sequence.nucleotide(find(variants_in_window.variant_index_start(j,1)==window_sequence.index));']);
                            eval(['window_sequence_var',num2str(n),'.variant_index_start(row_w,1) = variants_in_window.variant_index_start(j,1);']);
                            eval(['window_sequence_var',num2str(n),'.variant_index_end(row_w,1) = variants_in_window.variant_index_end(j,1);']);
                            eval(['window_sequence_var',num2str(n),'.variant_coord_start(row_w,1) = variants_in_window.variant_coord_start(j,1);']);
                            eval(['window_sequence_var',num2str(n),'.variant_coord_end(row_w,1) = variants_in_window.variant_coord_end(j,1);']);
                            eval(['window_sequence_var',num2str(n),'.rs_id(row_w,1) = variants_in_window.rs_id(j,1);']);
                            eval(['window_sequence_var',num2str(n),'.identifier(row_w,1) = variants_in_window.identifier(j,1);']);
                            eval(['window_sequence_var',num2str(n),'.variant_len(row_w,1) = variants_in_window.variant_len(j,1);']);
                            eval(['window_sequence_var',num2str(n),'.type{row_w,1} = ''Substitution'';']);
                            eval(['window_sequence_var',num2str(n),'.num_var{row_w,1} = ''Single Variant'';']);
                            row_w = row_w + 1;
                        elseif variants_in_window.identifier(j,1) == 2 & strcmp(strjoin(transpose(window_sequence.nucleotide(find(variants_in_window.variant_index_start(j,1)==window_sequence.index):find(variants_in_window.variant_index_end(j,1)==window_sequence.index))),''),variants_in_window.ref_snp(j,1)) == 1;
                            eval(['window_sequence_var',num2str(n),'.ref_snp(row_w,1) = variants_in_window.ref_snp(j,1);']);
                            eval(['window_sequence_var',num2str(n),'.alt_snp(row_w,1) = variants_in_window.alt_snp(j,1);']);
                            eval(['window_sequence_var',num2str(n),'.from_seq(row_w,1) = window_sequence.nucleotide(find(variants_in_window.variant_index_start(j,1)==window_sequence.index));']);
                            eval(['window_sequence_var',num2str(n),'.variant_index_start(row_w,1) = variants_in_window.variant_index_start(j,1);']);
                            eval(['window_sequence_var',num2str(n),'.variant_index_end(row_w,1) = variants_in_window.variant_index_end(j,1);']);
                            eval(['window_sequence_var',num2str(n),'.variant_coord_start(row_w,1) = variants_in_window.variant_coord_start(j,1);']);
                            eval(['window_sequence_var',num2str(n),'.variant_coord_end(row_w,1) = variants_in_window.variant_coord_end(j,1);']);
                            eval(['window_sequence_var',num2str(n),'.rs_id(row_w,1) = variants_in_window.rs_id(j,1);']);
                            eval(['window_sequence_var',num2str(n),'.identifier(row_w,1) = variants_in_window.identifier(j,1);']);
                            eval(['window_sequence_var',num2str(n),'.variant_len(row_w,1) = variants_in_window.variant_len(j,1);']);
                            eval(['window_sequence_var',num2str(n),'.type{row_w,1} = ''Deletion'';']);
                            eval(['window_sequence_var',num2str(n),'.num_var{row_w,1} = ''Single Variant'';']);
                            window_sequence_mod(find(variants_in_window.variant_index_start(j,1)==window_sequence.index):find(variants_in_window.variant_index_end(j,1)==window_sequence.index)) = {'D'};
                            ident_del_sites = strfind(window_sequence_mod,'D');
                            ident_del_sites = find(~cellfun(@isempty,ident_del_sites));
                            window_sequence_mod(ident_del_sites) = [];
                            eval(['window_sequence_var',num2str(n),'.window{row_w,1} = strjoin(transpose(window_sequence_mod),'''');']);
                            row_w = row_w + 1;
                        elseif variants_in_window.identifier(j,1) == 3
                            eval(['window_sequence_var',num2str(n),'.ref_snp(row_w,1) = variants_in_window.ref_snp(j,1);']);
                            eval(['window_sequence_var',num2str(n),'.alt_snp(row_w,1) = variants_in_window.alt_snp(j,1);']);
                            eval(['window_sequence_var',num2str(n),'.from_seq(row_w,1) = window_sequence.nucleotide(find(variants_in_window.variant_index_start(j,1)==window_sequence.index));']);
                            eval(['window_sequence_var',num2str(n),'.variant_index_start(row_w,1) = variants_in_window.variant_index_start(j,1);']);
                            eval(['window_sequence_var',num2str(n),'.variant_index_end(row_w,1) = variants_in_window.variant_index_end(j,1);']);
                            eval(['window_sequence_var',num2str(n),'.variant_coord_start(row_w,1) = variants_in_window.variant_coord_start(j,1);']);
                            eval(['window_sequence_var',num2str(n),'.variant_coord_end(row_w,1) = variants_in_window.variant_coord_end(j,1);']);
                            eval(['window_sequence_var',num2str(n),'.rs_id(row_w,1) = variants_in_window.rs_id(j,1);']);
                            eval(['window_sequence_var',num2str(n),'.identifier(row_w,1) = variants_in_window.identifier(j,1);']);
                            eval(['window_sequence_var',num2str(n),'.variant_len(row_w,1) = variants_in_window.variant_len(j,1);']);
                            eval(['window_sequence_var',num2str(n),'.type{row_w,1} = ''Insertion'';']);
                            eval(['window_sequence_var',num2str(n),'.num_var{row_w,1} = ''Single Variant'';']);
                            ins_length_increase = 0;
                            split_ins = transpose(cellstr(char(variants_in_window.alt_snp(j,1))')');
                            ins_length = length(split_ins);
                            index_ins_original = find(variants_in_window.variant_index_start(j,1)==window_sequence.index);
                            index_ins = find(variants_in_window.variant_index_start(j,1)==window_sequence.index)+ins_length_increase;
                            window_sequence_mod_final = window_sequence_mod(1:(index_ins-1));
                            window_sequence_mod_final((index_ins+ins_length):(length(window_sequence_mod)+ins_length)) = window_sequence_mod(index_ins:length(window_sequence_ins));
                            ins_length_increase = ins_length_increase + 1;
                            index_replace = index_ins;
                            ins_count = 0;
                            for i = 1:ins_length
                                window_sequence_mod_final(index_replace+ins_count) = split_ins(i,1);
                                ins_count = ins_count + 1;
                            end
                            eval(['window_sequence_var',num2str(n),'.window{row_w,1} = strjoin(transpose(window_sequence_mod_final),'''');']);
                            row_w = row_w + 1;
                        end
                    end
                end
            end
            
            if perform_variant_analysis == 1 && perform_wgs_analysis == 0
                row_var_index1 = 1;
                near_variant_index1 = [];
                if num_variants_window > 0
                    for w = 1:length(variants_in_window.variant_coord_start)
                        if isempty(find(abs(variants_in_window.variant_coord_start(w,1)-variants_in_window.variant_coord_start)<=(sgRNA_len+max(variants_in_window.variant_len)) & abs(variants_in_window.variant_coord_start(w,1)-variants_in_window.variant_coord_start) > 0)) == 0;
                            near_variant_index2 = find(abs(variants_in_window.variant_coord_start(w,1)-variants_in_window.variant_coord_start)<=(sgRNA_len+max(variants_in_window.variant_len)) & abs(variants_in_window.variant_coord_start(w,1)-variants_in_window.variant_coord_start) > 0);
                            near_variant_index1 = [near_variant_index1;near_variant_index2];
                            near_variant_index = unique(near_variant_index1);
                            find(abs(variants_in_window.variant_coord_start(w,1)-variants_in_window.variant_coord_start));
                            for z = 1:length(near_variant_index)
                                near_variants.ref_snp(z,1) = variants_in_window.ref_snp(near_variant_index(z,1));
                                near_variants.alt_snp(z,1) = variants_in_window.alt_snp(near_variant_index(z,1));
                                near_variants.variant_index_start(z,1) = variants_in_window.variant_index_start(near_variant_index(z,1));
                                near_variants.variant_index_end(z,1) = variants_in_window.variant_index_end(near_variant_index(z,1));
                                near_variants.variant_coord_start(z,1) = variants_in_window.variant_coord_start(near_variant_index(z,1));
                                near_variants.variant_coord_end(z,1) = variants_in_window.variant_coord_end(near_variant_index(z,1));
                                near_variants.rs_id(z,1) = variants_in_window.rs_id(near_variant_index(z,1));
                                near_variants.identifier(z,1) = variants_in_window.identifier(near_variant_index(z,1));
                                near_variants.variant_len(z,1) = variants_in_window.variant_len(near_variant_index(z,1));
                            end
                        end
                    end
                    if exist('near_variant_index','var') == 1
                        [near_variants.identifier type_index] = sort(variants_in_window.identifier,1);
                        near_variants.ref_snp = variants_in_window.ref_snp(type_index);
                        near_variants.alt_snp = variants_in_window.alt_snp(type_index);
                        near_variants.variant_index_start = variants_in_window.variant_index_start(type_index);
                        near_variants.variant_index_end = variants_in_window.variant_index_end(type_index);
                        near_variants.variant_coord_start = variants_in_window.variant_coord_start(type_index);
                        near_variants.variant_coord_end = variants_in_window.variant_coord_end(type_index);
                        near_variants.rs_id = variants_in_window.rs_id(type_index);
                        near_variants.variant_len = variants_in_window.variant_len(type_index);
                        num_ins_window = length(find(near_variants.identifier==3));
                        if length(near_variant_index) == 3 & exist('near_variants','var') == 1
                            num_iterations = length(near_variant_index)+1;
                            mult_var_guides_index1 = [1;3];
                            mult_var_guides_index2 = [1;2];
                            mult_var_guides_index3 = [2;3];
                            mult_var_guides_index = horzcat(mult_var_guides_index1,mult_var_guides_index2,mult_var_guides_index3);
                        elseif length(near_variant_index) == 2 & exist('near_variants','var') == 1
                            num_iterations = 1;
                        end
                        if length(near_variant_index) == 2 & exist('near_variants','var') == 1 | length(near_variant_index) == 3 & exist('near_variants','var') == 1
                            near_variants_holder = near_variants;
                            for y = 1:num_iterations
                                window_sequence_mult = window_sequence.nucleotide(:,1);
                                window_sequence_mult_index = window_sequence.index(:,1);
                                if length(near_variant_index) == 3 & y > 1
                                    clearvars near_variants
                                    near_variants.ref_snp = near_variants_holder.ref_snp(mult_var_guides_index(:,(y-1)));
                                    near_variants.alt_snp = near_variants_holder.alt_snp(mult_var_guides_index(:,(y-1)));
                                    near_variants.variant_index_start = near_variants_holder.variant_index_start(mult_var_guides_index(:,(y-1)));
                                    near_variants.variant_index_end = near_variants_holder.variant_index_end(mult_var_guides_index(:,(y-1)));
                                    near_variants.variant_coord_start = near_variants_holder.variant_coord_start(mult_var_guides_index(:,(y-1)));
                                    near_variants.variant_coord_end = near_variants_holder.variant_coord_end(mult_var_guides_index(:,(y-1)));
                                    near_variants.rs_id = near_variants_holder.rs_id(mult_var_guides_index(:,(y-1)));
                                    near_variants.identifier = near_variants_holder.identifier(mult_var_guides_index(:,(y-1)));
                                    near_variants.variant_len = near_variants_holder.variant_len(mult_var_guides_index(:,(y-1)));
                                end
                                if length(near_variant_index) == 3 & num_iterations > 1
                                    num_var_iterate = length(near_variants.ref_snp);
                                elseif num_iterations == 1;
                                    num_var_iterate = length(near_variant_index);
                                end
                                for h = 1:num_var_iterate
                                    if near_variants.identifier(h,1) == 1 & strcmp(window_sequence_mult(find(near_variants.variant_index_start(h,1)==window_sequence_mult_index)),near_variants.ref_snp(h,1)) == 1;
                                        window_sequence_mult(find(near_variants.variant_index_start(h,1)==window_sequence_mult_index),1) = near_variants.alt_snp(h,1);
                                        if isa(eval(['window_sequence_var',num2str(n),'.variant_index_start']),'double') == 1
                                            eval(['window_sequence_var',num2str(n),'.variant_index_start = num2cell(window_sequence_var',num2str(n),'.variant_index_start);']);
                                            eval(['window_sequence_var',num2str(n),'.variant_index_end = num2cell(window_sequence_var',num2str(n),'.variant_index_end);']);
                                            eval(['window_sequence_var',num2str(n),'.variant_coord_start = num2cell(window_sequence_var',num2str(n),'.variant_coord_start);']);
                                            eval(['window_sequence_var',num2str(n),'.variant_coord_end = num2cell(window_sequence_var',num2str(n),'.variant_coord_end);']);
                                            eval(['window_sequence_var',num2str(n),'.identifier = num2cell(window_sequence_var',num2str(n),'.identifier);']);
                                            eval(['window_sequence_var',num2str(n),'.variant_len = num2cell(window_sequence_var',num2str(n),'.variant_len);']);
                                        end
                                        if h == 1;
                                            eval(['window_sequence_var',num2str(n),'.ref_snp(row_w,1) = near_variants.ref_snp(h,1);']);
                                            eval(['window_sequence_var',num2str(n),'.alt_snp(row_w,1) = near_variants.alt_snp(h,1);']);
                                            eval(['window_sequence_var',num2str(n),'.from_seq(row_w,1) = window_sequence.nucleotide(find(near_variants.variant_index_start(h,1)==window_sequence.index));']);
                                            eval(['window_sequence_var',num2str(n),'.variant_index_start{row_w,1} = num2str(near_variants.variant_index_start(h,1));']);
                                            eval(['window_sequence_var',num2str(n),'.variant_index_end{row_w,1} = num2str(near_variants.variant_index_end(h,1));']);
                                            eval(['window_sequence_var',num2str(n),'.variant_coord_start{row_w,1} = num2str(near_variants.variant_coord_start(h,1));']);
                                            eval(['window_sequence_var',num2str(n),'.variant_coord_end{row_w,1} = num2str(near_variants.variant_coord_end(h,1));']);
                                            eval(['window_sequence_var',num2str(n),'.rs_id(row_w,1) = near_variants.rs_id(h,1);']);
                                            eval(['window_sequence_var',num2str(n),'.identifier{row_w,1} = num2str(near_variants.identifier(h,1));']);
                                            eval(['window_sequence_var',num2str(n),'.variant_len{row_w,1} = num2str(near_variants.variant_len(h,1));']);
                                            eval(['window_sequence_var',num2str(n),'.type{row_w,1} = ''Substitution'';']);
                                            eval(['window_sequence_var',num2str(n),'.num_var{row_w,1} = ''Multiple Variants'';']);
                                        end
                                        if h > 1
                                            eval(['window_sequence_var',num2str(n),'.ref_snp(row_w,1) = strcat(window_sequence_var',num2str(n),'.ref_snp(row_w,1),semicolon,near_variants.ref_snp(h,1));']);
                                            eval(['window_sequence_var',num2str(n),'.alt_snp(row_w,1) = strcat(window_sequence_var',num2str(n),'.alt_snp(row_w,1),semicolon,near_variants.alt_snp(h,1));']);
                                            eval(['window_sequence_var',num2str(n),'.from_seq(row_w,1) = strcat(window_sequence_var',num2str(n),'.from_seq(row_w,1),semicolon,window_sequence.nucleotide(find(near_variants.variant_index_start(h,1)==window_sequence.index)));']);
                                            eval(['window_sequence_var',num2str(n),'.type(row_w,1) = strcat(window_sequence_var',num2str(n),'.type(row_w,1),semicolon,''Substitution'');']);
                                            eval(['window_sequence_var',num2str(n),'.num_var(row_w,1) = strcat(window_sequence_var',num2str(n),'.num_var(row_w,1),semicolon,''Multiple Variants'');']);
                                            eval(['window_sequence_var',num2str(n),'.rs_id(row_w,1) = strcat(window_sequence_var',num2str(n),'.rs_id(row_w,1),semicolon,near_variants.rs_id(h,1));']);
                                            eval(['window_sequence_var',num2str(n),'.variant_index_start{row_w,1} = strcat(num2str(window_sequence_var',num2str(n),'.variant_index_start{row_w,1}),semicolon,num2str(near_variants.variant_index_start(h,1)));'])
                                            eval(['window_sequence_var',num2str(n),'.variant_index_end{row_w,1} = strcat(num2str(window_sequence_var',num2str(n),'.variant_index_end{row_w,1}),semicolon,num2str(near_variants.variant_index_end(h,1)));'])
                                            eval(['window_sequence_var',num2str(n),'.variant_coord_start{row_w,1} = strcat(num2str(window_sequence_var',num2str(n),'.variant_coord_start{row_w,1}),semicolon,num2str(near_variants.variant_coord_start(h,1)));'])
                                            eval(['window_sequence_var',num2str(n),'.variant_coord_end{row_w,1} = strcat(num2str(window_sequence_var',num2str(n),'.variant_coord_end{row_w,1}),semicolon,num2str(near_variants.variant_coord_end(h,1)));'])
                                            eval(['window_sequence_var',num2str(n),'.identifier{row_w,1} = strcat(num2str(window_sequence_var',num2str(n),'.identifier{row_w,1}),semicolon,num2str(near_variants.identifier(h,1)));'])
                                            eval(['window_sequence_var',num2str(n),'.variant_len{row_w,1} = strcat(num2str(window_sequence_var',num2str(n),'.variant_len{row_w,1}),semicolon,num2str(near_variants.variant_len(h,1)));'])
                                        end
                                    elseif near_variants.identifier(h,1) == 2 & strcmp(strjoin(transpose(window_sequence_mult(find(near_variants.variant_index_start(h,1)==window_sequence_mult_index):find(near_variants.variant_index_end(h,1)==window_sequence_mult_index))),''),near_variants.ref_snp(h,1)) == 1;
                                        if isa(eval(['window_sequence_var',num2str(n),'.variant_index_start']),'double') == 1
                                            eval(['window_sequence_var',num2str(n),'.variant_index_start = num2cell(window_sequence_var',num2str(n),'.variant_index_start);']);
                                            eval(['window_sequence_var',num2str(n),'.variant_index_end = num2cell(window_sequence_var',num2str(n),'.variant_index_end);']);
                                            eval(['window_sequence_var',num2str(n),'.variant_coord_start = num2cell(window_sequence_var',num2str(n),'.variant_coord_start);']);
                                            eval(['window_sequence_var',num2str(n),'.variant_coord_end = num2cell(window_sequence_var',num2str(n),'.variant_coord_end);']);
                                            eval(['window_sequence_var',num2str(n),'.identifier = num2cell(window_sequence_var',num2str(n),'.identifier);']);
                                            eval(['window_sequence_var',num2str(n),'.variant_len = num2cell(window_sequence_var',num2str(n),'.variant_len);']);
                                        end
                                        if h == 1;
                                            eval(['window_sequence_var',num2str(n),'.ref_snp(row_w,1) = near_variants.ref_snp(h,1);']);
                                            eval(['window_sequence_var',num2str(n),'.alt_snp(row_w,1) = near_variants.alt_snp(h,1);']);
                                            eval(['window_sequence_var',num2str(n),'.from_seq(row_w,1) = window_sequence.nucleotide(find(near_variants.variant_index_start(h,1)==window_sequence.index));']);
                                            eval(['window_sequence_var',num2str(n),'.variant_index_start{row_w,1} = num2str(near_variants.variant_index_start(h,1));']);
                                            eval(['window_sequence_var',num2str(n),'.variant_index_end{row_w,1} = num2str(near_variants.variant_index_end(h,1));']);
                                            eval(['window_sequence_var',num2str(n),'.variant_coord_start{row_w,1} = num2str(near_variants.variant_coord_start(h,1));']);
                                            eval(['window_sequence_var',num2str(n),'.variant_coord_end{row_w,1} = num2str(near_variants.variant_coord_end(h,1));']);
                                            eval(['window_sequence_var',num2str(n),'.rs_id(row_w,1) = near_variants.rs_id(h,1);']);
                                            eval(['window_sequence_var',num2str(n),'.identifier{row_w,1} = num2str(near_variants.identifier(h,1));']);
                                            eval(['window_sequence_var',num2str(n),'.variant_len{row_w,1} = num2str(near_variants.variant_len(h,1));']);
                                            eval(['window_sequence_var',num2str(n),'.type{row_w,1} = ''Deletion'';']);
                                            eval(['window_sequence_var',num2str(n),'.num_var{row_w,1} = ''Multiple Variants'';']);
                                        end
                                        if h > 1
                                            eval(['window_sequence_var',num2str(n),'.ref_snp(row_w,1) = strcat(window_sequence_var',num2str(n),'.ref_snp(row_w,1),semicolon,near_variants.ref_snp(h,1));']);
                                            eval(['window_sequence_var',num2str(n),'.alt_snp(row_w,1) = strcat(window_sequence_var',num2str(n),'.alt_snp(row_w,1),semicolon,near_variants.alt_snp(h,1));']);
                                            eval(['window_sequence_var',num2str(n),'.from_seq(row_w,1) = strcat(window_sequence_var',num2str(n),'.from_seq(row_w,1),semicolon,window_sequence.nucleotide(find(near_variants.variant_index_start(h,1)==window_sequence.index)));']);
                                            eval(['window_sequence_var',num2str(n),'.type(row_w,1) = strcat(window_sequence_var',num2str(n),'.type(row_w,1),semicolon,''Deletion'');']);
                                            eval(['window_sequence_var',num2str(n),'.num_var(row_w,1) = strcat(window_sequence_var',num2str(n),'.num_var(row_w,1),semicolon,''Multiple Variants'');']);
                                            eval(['window_sequence_var',num2str(n),'.rs_id(row_w,1) = strcat(window_sequence_var',num2str(n),'.rs_id(row_w,1),semicolon,near_variants.rs_id(h,1));']);
                                            eval(['window_sequence_var',num2str(n),'.variant_index_start{row_w,1} = strcat(num2str(window_sequence_var',num2str(n),'.variant_index_start{row_w,1}),semicolon,num2str(near_variants.variant_index_start(h,1)));'])
                                            eval(['window_sequence_var',num2str(n),'.variant_index_end{row_w,1} = strcat(num2str(window_sequence_var',num2str(n),'.variant_index_end{row_w,1}),semicolon,num2str(near_variants.variant_index_end(h,1)));'])
                                            eval(['window_sequence_var',num2str(n),'.variant_coord_start{row_w,1} = strcat(num2str(window_sequence_var',num2str(n),'.variant_coord_start{row_w,1}),semicolon,num2str(near_variants.variant_coord_start(h,1)));'])
                                            eval(['window_sequence_var',num2str(n),'.variant_coord_end{row_w,1} = strcat(num2str(window_sequence_var',num2str(n),'.variant_coord_end{row_w,1}),semicolon,num2str(near_variants.variant_coord_end(h,1)));'])
                                            eval(['window_sequence_var',num2str(n),'.identifier{row_w,1} = strcat(num2str(window_sequence_var',num2str(n),'.identifier{row_w,1}),semicolon,num2str(near_variants.identifier(h,1)));'])
                                            eval(['window_sequence_var',num2str(n),'.variant_len{row_w,1} = strcat(num2str(window_sequence_var',num2str(n),'.variant_len{row_w,1}),semicolon,num2str(near_variants.variant_len(h,1)));'])
                                        end
                                        window_sequence_mult(find(near_variants.variant_index_start(h,1)==window_sequence_mult_index):find(near_variants.variant_index_end(h,1)==window_sequence_mult_index)) = {'D'};
                                    elseif near_variants.identifier(h,1) == 3
                                        if isa(eval(['window_sequence_var',num2str(n),'.variant_index_start']),'double') == 1
                                            eval(['window_sequence_var',num2str(n),'.variant_index_start = num2cell(window_sequence_var',num2str(n),'.variant_index_start);']);
                                            eval(['window_sequence_var',num2str(n),'.variant_index_end = num2cell(window_sequence_var',num2str(n),'.variant_index_end);']);
                                            eval(['window_sequence_var',num2str(n),'.variant_coord_start = num2cell(window_sequence_var',num2str(n),'.variant_coord_start);']);
                                            eval(['window_sequence_var',num2str(n),'.variant_coord_end = num2cell(window_sequence_var',num2str(n),'.variant_coord_end);']);
                                            eval(['window_sequence_var',num2str(n),'.identifier = num2cell(window_sequence_var',num2str(n),'.identifier);']);
                                            eval(['window_sequence_var',num2str(n),'.variant_len = num2cell(window_sequence_var',num2str(n),'.variant_len);']);
                                        end
                                        if h == 1;
                                            eval(['window_sequence_var',num2str(n),'.ref_snp(row_w,1) = near_variants.ref_snp(h,1);']);
                                            eval(['window_sequence_var',num2str(n),'.alt_snp(row_w,1) = near_variants.alt_snp(h,1);']);
                                            eval(['window_sequence_var',num2str(n),'.from_seq(row_w,1) = window_sequence.nucleotide(find(near_variants.variant_index_start(h,1)==window_sequence.index));']);
                                            eval(['window_sequence_var',num2str(n),'.variant_index_start{row_w,1} = num2str(near_variants.variant_index_start(h,1));']);
                                            eval(['window_sequence_var',num2str(n),'.variant_index_end{row_w,1} = num2str(near_variants.variant_index_end(h,1));']);
                                            eval(['window_sequence_var',num2str(n),'.variant_coord_start{row_w,1} = num2str(near_variants.variant_coord_start(h,1));']);
                                            eval(['window_sequence_var',num2str(n),'.variant_coord_end{row_w,1} = num2str(near_variants.variant_coord_end(h,1));']);
                                            eval(['window_sequence_var',num2str(n),'.rs_id(row_w,1) = near_variants.rs_id(h,1);']);
                                            eval(['window_sequence_var',num2str(n),'.identifier{row_w,1} = num2str(near_variants.identifier(h,1));']);
                                            eval(['window_sequence_var',num2str(n),'.variant_len{row_w,1} = num2str(near_variants.variant_len(h,1));']);
                                            eval(['window_sequence_var',num2str(n),'.type{row_w,1} = ''Insertion'';']);
                                            eval(['window_sequence_var',num2str(n),'.num_var{row_w,1} = ''Multiple Variants'';']);
                                        end
                                        if h > 1
                                            eval(['window_sequence_var',num2str(n),'.ref_snp(row_w,1) = strcat(window_sequence_var',num2str(n),'.ref_snp(row_w,1),semicolon,near_variants.ref_snp(h,1));']);
                                            eval(['window_sequence_var',num2str(n),'.alt_snp(row_w,1) = strcat(window_sequence_var',num2str(n),'.alt_snp(row_w,1),semicolon,near_variants.alt_snp(h,1));']);
                                            eval(['window_sequence_var',num2str(n),'.from_seq(row_w,1) = strcat(window_sequence_var',num2str(n),'.from_seq(row_w,1),semicolon,window_sequence.nucleotide(find(near_variants.variant_index_start(h,1)==window_sequence.index)));']);
                                            eval(['window_sequence_var',num2str(n),'.type(row_w,1) = strcat(window_sequence_var',num2str(n),'.type(row_w,1),semicolon,''Insertion'');']);
                                            eval(['window_sequence_var',num2str(n),'.num_var(row_w,1) = strcat(window_sequence_var',num2str(n),'.num_var(row_w,1),semicolon,''Multiple Variants'');']);
                                            eval(['window_sequence_var',num2str(n),'.rs_id(row_w,1) = strcat(window_sequence_var',num2str(n),'.rs_id(row_w,1),semicolon,near_variants.rs_id(h,1));']);
                                            eval(['window_sequence_var',num2str(n),'.variant_index_start{row_w,1} = strcat(num2str(window_sequence_var',num2str(n),'.variant_index_start{row_w,1}),semicolon,num2str(near_variants.variant_index_start(h,1)));'])
                                            eval(['window_sequence_var',num2str(n),'.variant_index_end{row_w,1} = strcat(num2str(window_sequence_var',num2str(n),'.variant_index_end{row_w,1}),semicolon,num2str(near_variants.variant_index_end(h,1)));'])
                                            eval(['window_sequence_var',num2str(n),'.variant_coord_start{row_w,1} = strcat(num2str(window_sequence_var',num2str(n),'.variant_coord_start{row_w,1}),semicolon,num2str(near_variants.variant_coord_start(h,1)));'])
                                            eval(['window_sequence_var',num2str(n),'.variant_coord_end{row_w,1} = strcat(num2str(window_sequence_var',num2str(n),'.variant_coord_end{row_w,1}),semicolon,num2str(near_variants.variant_coord_end(h,1)));'])
                                            eval(['window_sequence_var',num2str(n),'.identifier{row_w,1} = strcat(num2str(window_sequence_var',num2str(n),'.identifier{row_w,1}),semicolon,num2str(near_variants.identifier(h,1)));'])
                                            eval(['window_sequence_var',num2str(n),'.variant_len{row_w,1} = strcat(num2str(window_sequence_var',num2str(n),'.variant_len{row_w,1}),semicolon,num2str(near_variants.variant_len(h,1)));'])
                                            
                                        end
                                        ins_length_increase = 0;
                                        index_ins_original = find(near_variants.variant_index_start(h,1)==window_sequence_mult_index);
                                        index_ins = find(near_variants.variant_index_start(h,1)==window_sequence_mult_index)+ins_length_increase;
                                        split_ins = transpose(cellstr(char(near_variants.alt_snp(h,1))')');
                                        ins_length = length(split_ins);
                                        window_sequence_mult1 = window_sequence_mult(1:(index_ins-1));
                                        window_sequence_mult1((index_ins+ins_length):(length(window_sequence_mult)+ins_length)) = window_sequence_mult(index_ins:length(window_sequence_mult));
                                        window_sequence_mult = window_sequence_mult1;
                                        window_sequence_mult_index1 = window_sequence_mult_index(1:(index_ins-1));
                                        window_sequence_mult_index1((index_ins+ins_length):(length(window_sequence_mult_index)+ins_length)) = window_sequence_mult_index(index_ins:length(window_sequence_mult_index));
                                        window_sequence_mult_index = window_sequence_mult_index1;
                                        ins_length_increase = ins_length_increase + 1;
                                        index_replace = index_ins;
                                        ins_count = 0;
                                        for i = 1:ins_length
                                            window_sequence_mult(index_replace+ins_count) = split_ins(i,1);
                                            ins_count = ins_count + 1;
                                        end
                                    end
                                end
                                if length(type_index) == 2 & exist('near_variants','var') == 1 | length(type_index) == 3 & exist('near_variants','var') == 1
                                    ident_del_sites = strfind(window_sequence_mult,'D');
                                    ident_del_sites = find(~cellfun(@isempty,ident_del_sites));
                                    window_sequence_mult(ident_del_sites) = [];
                                    eval(['window_sequence_var',num2str(n),'.window{row_w,1} = strjoin(transpose(window_sequence_mult),'''');']);
                                    row_w = row_w + 1;
                                end
                            end
                        end
                        clearvars sort_index index_ins near_variants_holder near_variant_index2 near_variant_index1 near_variant_index mult_var_guides_index
                    elseif exist('near_variant_index','var') == 0
                        eval(['window_sequence_var',num2str(n),'.variant_index_start = cellstr(num2str(window_sequence_var',num2str(n),'.variant_index_start));']);
                        eval(['window_sequence_var',num2str(n),'.variant_index_end = cellstr(num2str(window_sequence_var',num2str(n),'.variant_index_end));']);
                        eval(['window_sequence_var',num2str(n),'.variant_coord_start = cellstr(num2str(window_sequence_var',num2str(n),'.variant_coord_start));']);
                        eval(['window_sequence_var',num2str(n),'.variant_coord_end = cellstr(num2str(window_sequence_var',num2str(n),'.variant_coord_end));']);
                        eval(['window_sequence_var',num2str(n),'.identifier = cellstr(num2str(window_sequence_var',num2str(n),'.identifier));']);
                        eval(['window_sequence_var',num2str(n),'.variant_len = cellstr(num2str(window_sequence_var',num2str(n),'.variant_len));']);
                    end
                end
            end
        end
        
        orig_guide_index1 = [];
        orig_guide_index2 = [];
        orig_guide_index1_window_start1 = [];
        orig_guide_index1_window_end1 = [];
        orig_guide_index2_window_start1 = [];
        orig_guide_index2_window_end1 = [];
        for n = 1:num_windows
            eval(['original_window_seq = window_sequence_var',num2str(n),'.window(1,1);']);
            eval(['original_window_index_start = window_sequence_var',num2str(n),'.window_start;']);
            eval(['original_window_index_end = window_sequence_var',num2str(n),'.window_end;']);
            original_window_seq = char(original_window_seq);
            if length(PAM_sequence_cell)>1
                for u  = 1:length(PAM_sequence_cell)
                    PAM_sequence_cell_rc{u,1} = seqrcomplement(PAM_sequence_cell{u,1});
                end
            else PAM_sequence_cell_rc = cellstr(seqrcomplement(PAM_sequence_cell{:}));
            end
            orig_guide_index1 = [];
            orig_guide_index2 = [];
            for m = 1:length(PAM_sequence_cell)
                orig_guide_index1_temp_index = transpose(strfind(original_window_seq,PAM_sequence_cell{m,1}));
                orig_guide_index1_temp = orig_guide_index1_temp_index + original_window_index_start - 1;
                orig_guide_index1 = [orig_guide_index1;orig_guide_index1_temp];
            end
            if five_prime_PAM == 0
                orig_guide_index1_window_start1 = [orig_guide_index1_window_start1;(orig_guide_index1-(N_seq_in_PAM+sgRNA_len))];
                orig_guide_index1_window_end1 = [orig_guide_index1_window_end1;(orig_guide_index1+(PAM_len-N_seq_in_PAM-1))];
            end
            if five_prime_PAM == 1
                orig_guide_index1_window_start1 = [orig_guide_index1_window_start1;orig_guide_index1];
                orig_guide_index1_window_end1 = [orig_guide_index1_window_end1;(orig_guide_index1+(PAM_len+sgRNA_len-1))];
            end
            
            
            for l = 1:length(PAM_sequence_cell)
                orig_guide_index2_temp_index = transpose(strfind(original_window_seq,PAM_sequence_cell_rc{l,1}));
                orig_guide_index2_temp = orig_guide_index2_temp_index + original_window_index_start - 1;
                orig_guide_index2 = [orig_guide_index2;orig_guide_index2_temp];
            end
            if five_prime_PAM == 0
                orig_guide_index2_window_start1 = [orig_guide_index2_window_start1;orig_guide_index2];
                orig_guide_index2_window_end1 = [orig_guide_index2_window_end1;(orig_guide_index2+(sgRNA_len+PAM_len-1))];
            end
            if five_prime_PAM == 1
                orig_guide_index2_window_start1 = [orig_guide_index2_window_start1;orig_guide_index2-(N_seq_in_PAM+sgRNA_len)];
                orig_guide_index2_window_end1 = [orig_guide_index2_window_end1;(orig_guide_index2+PAM_len-N_seq_in_PAM-1)];
            end
            
            clearvars orig_guide_index1_temp_index orig_guide_index1_temp orig_guide_index1 original_window_seq original_window_seq_rc
            clearvars orig_guide_index2_temp_index orig_guide_index2_temp orig_guide_index2 original_window_seq original_window_seq_rc
        end
        
        row = 1;
        guide_ID_num = 1;
        
        if isempty(orig_guide_index1_window_start1) == 0
            for n = 1:length(orig_guide_index1_window_start1)
                zeros_to_add = [];
                if five_prime_PAM == 0 && orig_guide_index1_window_start1(n,1) > 0 && orig_guide_index1_window_end1(n,1)<=length(reference_seq)
                    orig_guide_index1_window_start(row,1) = orig_guide_index1_window_start1(n,1);
                    orig_guide_index1_window_end(row,1) = orig_guide_index1_window_end1(n,1);
                    orig_guide_index1_window_start_coord(row,1) = reference_sequence.coord(orig_guide_index1_window_start(row,1));
                    orig_guide_index1_window_end_coord(row,1) = reference_sequence.coord(orig_guide_index1_window_end(row,1));
                    orig_guide_index1_window_dsb_index(row,1) = orig_guide_index1_window_start(row,1)+(sgRNA_len-4);
                    orig_guide_index1_window_dsb_coord(row,1) = reference_sequence.coord(orig_guide_index1_window_dsb_index(row,1));
                    row = row + 1;
                elseif five_prime_PAM == 1 && orig_guide_index1_window_start1(n,1) > 0 && orig_guide_index1_window_end1(n,1)<=length(reference_seq)
                    orig_guide_index1_window_start(row,1) = orig_guide_index1_window_start1(n,1);
                    orig_guide_index1_window_end(row,1) = orig_guide_index1_window_end1(n,1);
                    orig_guide_index1_window_start_coord(row,1) = reference_sequence.coord(orig_guide_index1_window_start(row,1));
                    orig_guide_index1_window_end_coord(row,1) = reference_sequence.coord(orig_guide_index1_window_end(row,1));
                    orig_guide_index1_window_dsb_index(row,1) = orig_guide_index1_window_start(row,1)+(PAM_len+18);
                    orig_guide_index1_window_dsb_coord(row,1) = reference_sequence.coord(orig_guide_index1_window_dsb_index(row,1));
                    row = row + 1;
                end
            end
        end
        row = 1;
        if isempty(orig_guide_index2_window_start1) == 0
            for n = 1:length(orig_guide_index2_window_start1)
                if five_prime_PAM == 0 && orig_guide_index2_window_start1(n,1) > 0 && orig_guide_index2_window_end1(n,1)<=length(reference_seq)
                    orig_guide_index2_window_end(row,1) = orig_guide_index2_window_start1(n,1);
                    orig_guide_index2_window_start(row,1) = orig_guide_index2_window_end1(n,1);
                    orig_guide_index2_window_end_coord(row,1) = reference_sequence.coord(orig_guide_index2_window_end(row,1));
                    orig_guide_index2_window_start_coord(row,1) = reference_sequence.coord(orig_guide_index2_window_start(row,1));
                    orig_guide_index2_window_dsb_index(row,1) = orig_guide_index2_window_end(row,1)+(PAM_len+2);
                    orig_guide_index2_window_dsb_coord(row,1) = reference_sequence.coord(orig_guide_index2_window_dsb_index(row,1));
                    row = row + 1;
                elseif five_prime_PAM == 1 && orig_guide_index2_window_start1(n,1) > 0 && orig_guide_index2_window_end1(n,1)<=length(reference_seq)
                    orig_guide_index2_window_end(row,1) = orig_guide_index2_window_start1(n,1);
                    orig_guide_index2_window_start(row,1) = orig_guide_index2_window_end1(n,1);
                    orig_guide_index2_window_end_coord(row,1) = reference_sequence.coord(orig_guide_index2_window_end(row,1));
                    orig_guide_index2_window_start_coord(row,1) = reference_sequence.coord(orig_guide_index2_window_start(row,1));
                    orig_guide_index2_window_dsb_index(row,1) = orig_guide_index2_window_start(row,1)-(PAM_len+18);
                    orig_guide_index2_window_dsb_coord(row,1) = reference_sequence.coord(orig_guide_index2_window_dsb_index(row,1));
                    row = row + 1;
                end
                
            end
        end
        clearvars no_variant_analysis no_variant_analysis_rc
        
        if exist('orig_guide_index1_window_start','var') == 1 & five_prime_PAM == 0;
            for n = 1:length(orig_guide_index1_window_end)
                no_variant_analysis.guides{n,1} = strjoin(transpose(reference_sequence.nucleotide(orig_guide_index1_window_start(n,1):orig_guide_index1_window_end(n,1))),'');
                no_variant_analysis.guides_only{n,1} = no_variant_analysis.guides{n,1}(1:sgRNA_len);
                no_variant_analysis.pam_only{n,1} = no_variant_analysis.guides{n,1}(sgRNA_len+1:length(no_variant_analysis.guides{n,1}));
                no_variant_analysis.strand{n,1} = '+';
                no_variant_analysis.dsb_index(n,1) = orig_guide_index1_window_start(n,1)+(sgRNA_len-4);
                no_variant_analysis.dsb_coord(n,1) = reference_sequence.coord(no_variant_analysis.dsb_index(n,1));
                no_variant_analysis.guide_index_start(n,1) = orig_guide_index1_window_start(n,1);
                no_variant_analysis.guide_index_end(n,1) = orig_guide_index1_window_end(n,1);
                no_variant_analysis.guide_coord_start(n,1) = orig_guide_index1_window_start_coord(n,1);
                no_variant_analysis.guide_coord_end(n,1) = orig_guide_index1_window_end_coord(n,1);
            end
            clearvars orig_guide_index1_window_start orig_guide_index1_window_end orig_guide_index1_window_start_coord orig_guide_index1_window_end_coord
        end
        if exist('orig_guide_index2_window_end','var') == 1 & five_prime_PAM == 0;
            for n = 1:length(orig_guide_index2_window_end)
                no_variant_analysis_rc.guides_rc{n,1} = strjoin(transpose(reference_sequence.nucleotide(orig_guide_index2_window_end(n,1):orig_guide_index2_window_start(n,1))),'');
                no_variant_analysis_rc.guides{n,1} = seqrcomplement(no_variant_analysis_rc.guides_rc{n,1});
                no_variant_analysis_rc.guides_only{n,1} = no_variant_analysis_rc.guides{n,1}(1:sgRNA_len);
                no_variant_analysis_rc.pam_only{n,1} = no_variant_analysis_rc.guides{n,1}(sgRNA_len+1:length(no_variant_analysis_rc.guides{n,1}));
                no_variant_analysis_rc.strand{n,1} = '-';
                no_variant_analysis_rc.dsb_index(n,1) = orig_guide_index2_window_end(n,1)+(PAM_len+3);
                no_variant_analysis_rc.dsb_coord(n,1) = reference_sequence.coord(no_variant_analysis_rc.dsb_index(n,1));
                no_variant_analysis_rc.guide_index_start(n,1) = orig_guide_index2_window_start(n,1);
                no_variant_analysis_rc.guide_index_end(n,1)= orig_guide_index2_window_end(n,1);
                no_variant_analysis_rc.guide_coord_start(n,1) = orig_guide_index2_window_start_coord(n,1);
                no_variant_analysis_rc.guide_coord_end(n,1) = orig_guide_index2_window_end_coord(n,1);
            end
            clearvars orig_guide_index2_window_start orig_guide_index2_window_end orig_guide_index2_window_start_coord orig_guide_index2_window_end_coord
        end
        
        if exist('orig_guide_index1_window_start','var') == 1 & five_prime_PAM == 1;
            for n = 1:length(orig_guide_index1_window_end)
                no_variant_analysis.guides{n,1} = strjoin(transpose(reference_sequence.nucleotide(orig_guide_index1_window_start(n,1):orig_guide_index1_window_end(n,1))),'');
                no_variant_analysis.guides_only{n,1} = no_variant_analysis.guides{n,1}((PAM_len+1):(PAM_len+sgRNA_len));
                no_variant_analysis.pam_only{n,1} = no_variant_analysis.guides{n,1}(1:PAM_len);
                no_variant_analysis.strand{n,1} = '+';
                no_variant_analysis.dsb_index(n,1) = orig_guide_index1_window_start(n,1)+(PAM_len+18);
                no_variant_analysis.dsb_coord(n,1) = reference_sequence.coord(no_variant_analysis.dsb_index(n,1));
                no_variant_analysis.guide_index_start(n,1) = orig_guide_index1_window_start(n,1);
                no_variant_analysis.guide_index_end(n,1) = orig_guide_index1_window_end(n,1);
                no_variant_analysis.guide_coord_start(n,1) = orig_guide_index1_window_start_coord(n,1);
                no_variant_analysis.guide_coord_end(n,1) = orig_guide_index1_window_end_coord(n,1);
            end
            clearvars orig_guide_index1_window_start orig_guide_index1_window_end orig_guide_index1_window_start_coord orig_guide_index1_window_end_coord
        end
        if exist('orig_guide_index2_window_end','var') == 1 & five_prime_PAM == 1;
            for n = 1:length(orig_guide_index2_window_end)
                no_variant_analysis_rc.guides_rc{n,1} = strjoin(transpose(reference_sequence.nucleotide(orig_guide_index2_window_end(n,1):orig_guide_index2_window_start(n,1))),'');
                no_variant_analysis_rc.guides{n,1} = seqrcomplement(no_variant_analysis_rc.guides_rc{n,1});
                no_variant_analysis_rc.guides_only{n,1} = no_variant_analysis_rc.guides{n,1}((PAM_len+1):(PAM_len+sgRNA_len));
                no_variant_analysis_rc.pam_only{n,1} = no_variant_analysis_rc.guides{n,1}(1:PAM_len);
                no_variant_analysis_rc.strand{n,1} = '-';
                no_variant_analysis_rc.dsb_index(n,1) = orig_guide_index2_window_start(n,1)-(PAM_len+18);
                no_variant_analysis_rc.dsb_coord(n,1) = reference_sequence.coord(no_variant_analysis_rc.dsb_index(n,1));
                no_variant_analysis_rc.guide_index_start(n,1) = orig_guide_index2_window_start(n,1);
                no_variant_analysis_rc.guide_index_end(n,1)= orig_guide_index2_window_end(n,1);
                no_variant_analysis_rc.guide_coord_start(n,1) = orig_guide_index2_window_start_coord(n,1);
                no_variant_analysis_rc.guide_coord_end(n,1) = orig_guide_index2_window_end_coord(n,1);
            end
            clearvars orig_guide_index2_window_start orig_guide_index2_window_end orig_guide_index2_window_start_coord orig_guide_index2_window_end_coord
        end
        
        if exist('no_variant_analysis','var') == 1 & exist('no_variant_analysis_rc','var') == 1
            original_guides_full_temp.guides = [no_variant_analysis.guides;no_variant_analysis_rc.guides];
            original_guides_full_temp.guides_only = [no_variant_analysis.guides_only;no_variant_analysis_rc.guides_only];
            original_guides_full_temp.pam_only = [no_variant_analysis.pam_only;no_variant_analysis_rc.pam_only];
            original_guides_full_temp.strand = [no_variant_analysis.strand;no_variant_analysis_rc.strand];
            original_guides_full_temp.dsb_index = [no_variant_analysis.dsb_index;no_variant_analysis_rc.dsb_index];
            original_guides_full_temp.dsb_coord = [no_variant_analysis.dsb_coord;no_variant_analysis_rc.dsb_coord];
            original_guides_full_temp.guide_index_start = [no_variant_analysis.guide_index_start;no_variant_analysis_rc.guide_index_start];
            original_guides_full_temp.guide_index_end = [no_variant_analysis.guide_index_end;no_variant_analysis_rc.guide_index_end];
            original_guides_full_temp.guide_coord_start = [no_variant_analysis.guide_coord_start;no_variant_analysis_rc.guide_coord_start];
            original_guides_full_temp.guide_coord_end = [no_variant_analysis.guide_coord_end;no_variant_analysis_rc.guide_coord_end];
        elseif exist('no_variant_analysis','var') == 1 & exist('no_variant_analysis_rc','var') == 0
            original_guides_full_temp.guides = no_variant_analysis.guides;
            original_guides_full_temp.guides_only = no_variant_analysis.guides_only;
            original_guides_full_temp.pam_only = no_variant_analysis.pam_only;
            original_guides_full_temp.strand = no_variant_analysis.strand;
            original_guides_full_temp.dsb_index = no_variant_analysis.dsb_index;
            original_guides_full_temp.dsb_coord = no_variant_analysis.dsb_coord;
            original_guides_full_temp.guide_index_start = no_variant_analysis.guide_index_start;
            original_guides_full_temp.guide_index_end = no_variant_analysis.guide_index_end;
            original_guides_full_temp.guide_coord_start = no_variant_analysis.guide_coord_start;
            original_guides_full_temp.guide_coord_end = no_variant_analysis.guide_coord_end;
        elseif exist('no_variant_analysis','var') == 0 & exist('no_variant_analysis_rc','var') == 1
            original_guides_full_temp.guides = no_variant_analysis_rc.guides;
            original_guides_full_temp.guides_only = no_variant_analysis_rc.guides_only;
            original_guides_full_temp.pam_only = no_variant_analysis_rc.pam_only;
            original_guides_full_temp.strand = no_variant_analysis_rc.strand;
            original_guides_full_temp.dsb_index = no_variant_analysis_rc.dsb_index;
            original_guides_full_temp.dsb_coord = no_variant_analysis_rc.dsb_coord;
            original_guides_full_temp.guide_index_start = no_variant_analysis_rc.guide_index_start;
            original_guides_full_temp.guide_index_end = no_variant_analysis_rc.guide_index_end;
            original_guides_full_temp.guide_coord_start = no_variant_analysis_rc.guide_coord_start;
            original_guides_full_temp.guide_coord_end = no_variant_analysis_rc.guide_coord_end;
        end
        
        clearvars no_variant_analysis no_variant_analysis_rc
        if exist('original_guides_full_temp','var') == 1
            [original_guides_full.guides unique_index_full] = unique(original_guides_full_temp.guides);
            original_guides_full.guides_only = original_guides_full_temp.guides_only(unique_index_full);
            original_guides_full.pam_only = original_guides_full_temp.pam_only(unique_index_full);
            original_guides_full.strand = original_guides_full_temp.strand(unique_index_full);
            original_guides_full.dsb_index = original_guides_full_temp.dsb_index(unique_index_full);
            original_guides_full.dsb_coord = original_guides_full_temp.dsb_coord(unique_index_full);
            original_guides_full.guide_index_start = original_guides_full_temp.guide_index_start(unique_index_full);
            original_guides_full.guide_index_end = original_guides_full_temp.guide_index_end(unique_index_full);
            original_guides_full.guide_coord_start = original_guides_full_temp.guide_coord_start(unique_index_full);
            original_guides_full.guide_coord_end = original_guides_full_temp.guide_coord_end(unique_index_full);
            clearvars original_guides_full_temp
        end
        
        if twenty_extra_bases == 1 && perform_haplotype_analysis == 0 && exist('original_guides_full','var') == 1;
            edge_index_left = find(original_guides_full.dsb_coord<(seq_bpstart+sgRNA_len));
            original_guides_full.guides(edge_index_left) = [];
            original_guides_full.guides_only(edge_index_left) = [];
            original_guides_full.pam_only(edge_index_left) = [];
            original_guides_full.strand(edge_index_left) = [];
            original_guides_full.dsb_index(edge_index_left) = [];
            original_guides_full.dsb_coord(edge_index_left) = [];
            original_guides_full.guide_index_start(edge_index_left) = [];
            original_guides_full.guide_index_end(edge_index_left) = [];
            original_guides_full.guide_coord_start(edge_index_left) = [];
            original_guides_full.guide_coord_end(edge_index_left) = [];
            edge_index_right = find(original_guides_full.dsb_coord>(seq_bpend-sgRNA_len));
            original_guides_full.guides(edge_index_right) = [];
            original_guides_full.guides_only(edge_index_right) = [];
            original_guides_full.pam_only(edge_index_right) = [];
            original_guides_full.strand(edge_index_right) = [];
            original_guides_full.dsb_index(edge_index_right) = [];
            original_guides_full.dsb_coord(edge_index_right) = [];
            original_guides_full.guide_index_start(edge_index_right) = [];
            original_guides_full.guide_index_end(edge_index_right) = [];
            original_guides_full.guide_coord_start(edge_index_right) = [];
            original_guides_full.guide_coord_end(edge_index_right) = [];
        end
        
        if perform_haplotype_analysis == 0 && perform_variant_analysis == 0 && perform_wgs_analysis == 0 && exist('original_guides_full','var') == 1;
            if a == starting_number;
                for n = 1:length(original_guides_full.guides)
                    zeros_to_add = [];
                    if exist('guide_ID_num_temp','var') == 0
                        guide_ID_num_temp = n;
                    else guide_ID_num_temp = guide_ID_num_temp + 1;
                    end
                    num_zeros_add = length(num2str(guide_ID_num_temp));
                    if num_zeros_add < max_number_guide_ID_digits
                        for j = 1:(max_number_guide_ID_digits-num_zeros_add)
                            zeros_to_add = strcat(zeros_to_add,'0');
                        end
                        original_guides_full.guide_ID_number{n,1} = strcat('sgRNA_',zeros_to_add,num2str(guide_ID_num_temp));
                    elseif num_zeros_add == max_number_guide_ID_digits
                        original_guides_full.guide_ID_number{n,1} = ['sgRNA_',num2str(guide_ID_num_temp),''];
                    end
                end
            elseif a > starting_number
                for n = 1:length(original_guides_full.guides)
                    zeros_to_add = [];
                    guide_ID_num_temp = guide_ID_num_temp+1;
                    num_zeros_add = length(num2str(guide_ID_num_temp));
                    if num_zeros_add < max_number_guide_ID_digits
                        for j = 1:(max_number_guide_ID_digits-num_zeros_add)
                            zeros_to_add = strcat(zeros_to_add,'0');
                        end
                        original_guides_full.guide_ID_number{n,1} = strcat('sgRNA_',zeros_to_add,num2str(guide_ID_num_temp));
                    elseif num_zeros_add == max_number_guide_ID_digits
                        original_guides_full.guide_ID_number{n,1} = ['sgRNA_',num2str(guide_ID_num_temp),''];
                    end
                end
                
            end
        end
        
        %% Multiple Match Analysis
        
        if exist('original_guides_full','var') == 1 && multiple_match_analysis == 1
            for n = 1:length(original_guides_full.guides)
                temp_guide = seqrcomplement(original_guides_full.guides_only{n,1});
                if isempty(strfind(reference_seq,original_guides_full.guides_only{n,1})) == 0 && isempty(strfind(reference_seq,temp_guide)) == 0;
                    original_guides_full.mult_match_count(n,1) = size(strfind(reference_seq,original_guides_full.guides_only{n,1}),2) + size(strfind(reference_seq,temp_guide),2);
                    category = 1;
                elseif isempty(strfind(reference_seq,original_guides_full.guides_only{n,1})) == 0 && isempty(strfind(reference_seq,temp_guide)) == 1
                    original_guides_full.mult_match_count(n,1) = size(strfind(reference_seq,original_guides_full.guides{n,1}),2);
                    category = 2;
                elseif isempty(strfind(reference_seq,original_guides_full.guides_only{n,1})) == 1 && isempty(strfind(reference_seq,temp_guide)) == 0
                    original_guides_full.mult_match_count(n,1) = size(strfind(reference_seq,temp_guide),2);
                    category = 3;
                end
                if original_guides_full.mult_match_count(n,1) > 1
                    original_guides_full.category(n,1) = category;
                else original_guides_full.category(n,1) = 0;
                end
            end
            
            row = 1;
            mult_match_index = find(original_guides_full.mult_match_count>1);
            for n = 1:length(mult_match_index)
                if five_prime_PAM == 0;
                    if original_guides_full.category(mult_match_index(n,1)) == 1
                        num_duplicates = transpose(strfind(reference_seq,char(original_guides_full.guides_only(mult_match_index(n,1)))));
                        num_duplicates_rc = transpose(strfind(reference_seq,seqrcomplement(char(original_guides_full.guides_only(mult_match_index(n,1))))));
                        for i = 1:length(num_duplicates)
                            multiple_match.guides{row,1} = char(reference_seq(num_duplicates(i,1):(num_duplicates(i,1)+(sgRNA_len-1)+PAM_len)));
                            multiple_match.guides_only{row,1} = multiple_match.guides{row,1}(1:sgRNA_len);
                            multiple_match.pam_only{row,1} = multiple_match.guides{row,1}(sgRNA_len+1:length(multiple_match.guides{row,1}));
                            multiple_match.pam(row,1) = {PAM_sequence1};
                            multiple_match.strand{row,1} = '+';
                            multiple_match.dsb_index(row,1) = num_duplicates(i,1)+(sgRNA_len-4);
                            multiple_match.dsb_coord(row,1) = reference_sequence.coord(multiple_match.dsb_index(row,1));
                            multiple_match.guide_index_start(row,1) = num_duplicates(i,1);
                            multiple_match.guide_index_end(row,1) = num_duplicates(i,1)+(sgRNA_len-1)+PAM_len;
                            multiple_match.guide_coord_start(row,1) = reference_sequence.coord(num_duplicates(i,1));
                            multiple_match.guide_coord_end(row,1) = reference_sequence.coord(num_duplicates(i,1)+(sgRNA_len-1)+PAM_len);
                            row = row + 1;
                        end
                        for i = 1:length(num_duplicates_rc)
                            multiple_match.guides{row,1} = seqrcomplement(char(reference_seq((num_duplicates_rc(i,1)-PAM_len):(num_duplicates_rc(i,1)+(sgRNA_len-1)))));
                            multiple_match.guides_only{row,1} = multiple_match.guides{row,1}(1:sgRNA_len);
                            multiple_match.pam_only{row,1} = multiple_match.guides{row,1}(sgRNA_len+1:length(multiple_match.guides{row,1}));
                            multiple_match.pam(row,1) = {PAM_sequence1};
                            multiple_match.strand{row,1} = '-';
                            multiple_match.dsb_index(row,1) = num_duplicates_rc(i,1)+PAM_len+3;
                            multiple_match.dsb_coord(row,1) = reference_sequence.coord(multiple_match.dsb_index(row,1));
                            multiple_match.guide_index_end(row,1) = num_duplicates_rc(i,1);
                            multiple_match.guide_index_start(row,1) = num_duplicates_rc(i,1)+(sgRNA_len-1)+PAM_len;
                            multiple_match.guide_coord_end(row,1) = reference_sequence.coord(num_duplicates_rc(i,1));
                            multiple_match.guide_coord_start(row,1) = reference_sequence.coord(num_duplicates_rc(i,1)+(sgRNA_len-1)+PAM_len);
                            row = row + 1;
                        end
                    elseif original_guides_full.category(mult_match_index(n,1)) == 2
                        num_duplicates = transpose(strfind(reference_seq,char(original_guides_full.guides_only(mult_match_index(n,1)))));
                        for i = 1:length(num_duplicates)
                            multiple_match.guides{row,1} = char(reference_seq(num_duplicates(i,1):(num_duplicates(i,1)+(sgRNA_len-1)+PAM_len)));
                            multiple_match.guides_only{row,1} = multiple_match.guides{row,1}(1:sgRNA_len);
                            multiple_match.pam_only{row,1} = multiple_match.guides{row,1}(sgRNA_len+1:length(multiple_match.guides{row,1}));
                            multiple_match.pam(row,1) = {PAM_sequence1};
                            multiple_match.strand{row,1} = '+';
                            multiple_match.dsb_index(row,1) = num_duplicates(i,1)+(sgRNA_len-4);
                            multiple_match.dsb_coord(row,1) = reference_sequence.coord(multiple_match.dsb_index(row,1));
                            multiple_match.guide_index_start(row,1) = num_duplicates(i,1);
                            multiple_match.guide_index_end(row,1) = num_duplicates(i,1)+(sgRNA_len-1)+PAM_len;
                            multiple_match.guide_coord_start(row,1) = reference_sequence.coord(num_duplicates(i,1));
                            multiple_match.guide_coord_end(row,1) = reference_sequence.coord(num_duplicates(i,1)+(sgRNA_len-1)+PAM_len);
                            row = row + 1;
                        end
                    elseif original_guides_full.category(mult_match_index(n,1)) == 3
                        num_duplicates_rc = transpose(strfind(reference_seq,seqrcomplement(char(original_guides_full.guides_only(mult_match_index(n,1))))));
                        for i = 1:length(num_duplicates_rc)
                            multiple_match.guides{row,1} = seqrcomplement(char(reference_seq((num_duplicates_rc(i,1)-PAM_len):(num_duplicates_rc(i,1)+(sgRNA_len-1)))));
                            multiple_match.guides_only{row,1} = multiple_match.guides{row,1}(1:sgRNA_len);
                            multiple_match.pam_only{row,1} = multiple_match.guides{row,1}(sgRNA_len+1:length(multiple_match.guides{row,1}));
                            multiple_match.pam(row,1) = {PAM_sequence1};
                            multiple_match.strand{row,1} = '-';
                            multiple_match.dsb_index(row,1) = num_duplicates_rc(i,1)+PAM_len+3;
                            multiple_match.dsb_coord(row,1) = reference_sequence.coord(multiple_match.dsb_index(row,1));
                            multiple_match.guide_index_end(row,1) = num_duplicates_rc(i,1);
                            multiple_match.guide_index_start(row,1) = num_duplicates_rc(i,1)+(sgRNA_len-1)+PAM_len;
                            multiple_match.guide_coord_end(row,1) = reference_sequence.coord(num_duplicates_rc(i,1));
                            multiple_match.guide_coord_start(row,1) = reference_sequence.coord(num_duplicates_rc(i,1)+(sgRNA_len-1)+PAM_len);
                            row = row + 1;
                        end
                    end
                elseif five_prime_PAM == 1;
                    if original_guides_full.category(mult_match_index(n,1)) == 1
                        num_duplicates = transpose(strfind(reference_seq,char(original_guides_full.guides_only(mult_match_index(n,1)))));
                        num_duplicates_rc = transpose(strfind(reference_seq,seqrcomplement(char(original_guides_full.guides_only(mult_match_index(n,1))))));
                        for i = 1:length(num_duplicates)
                            multiple_match.guides{row,1} = char(reference_seq((num_duplicates(i,1)-PAM_len):(num_duplicates(i,1)+(sgRNA_len-1))));
                            multiple_match.guides_only{row,1} = multiple_match.guides{row,1}((PAM_len+1):(PAM_len+sgRNA_len));
                            multiple_match.pam_only{row,1} = multiple_match.guides{row,1}(1:PAM_len);
                            multiple_match.pam(row,1) = {PAM_sequence1};
                            multiple_match.strand{row,1} = '+';
                            multiple_match.dsb_index(row,1) = num_duplicates(i,1)+18;
                            multiple_match.dsb_coord(row,1) = reference_sequence.coord(multiple_match.dsb_index(row,1));
                            multiple_match.guide_index_start(row,1) = num_duplicates(i,1)-PAM_len;
                            multiple_match.guide_index_end(row,1) = num_duplicates(i,1)+(sgRNA_len-1);
                            multiple_match.guide_coord_start(row,1) = reference_sequence.coord(num_duplicates(i,1)-PAM_len);
                            multiple_match.guide_coord_end(row,1) = reference_sequence.coord(num_duplicates(i,1)+(sgRNA_len-1));
                            row = row + 1;
                        end
                        for i = 1:length(num_duplicates_rc)
                            multiple_match.guides{row,1} = seqrcomplement(char(reference_seq(num_duplicates_rc(i,1):(num_duplicates_rc(i,1)+(sgRNA_len-1)+PAM_len))));
                            multiple_match.guides_only{row,1} = multiple_match.guides{row,1}(PAM_len+1:length(multiple_match.guides{row,1}));
                            multiple_match.pam_only{row,1} = multiple_match.guides{row,1}(1:PAM_len);
                            multiple_match.pam(row,1) = {PAM_sequence1};
                            multiple_match.strand{row,1} = '-';
                            multiple_match.guide_index_end(row,1) = num_duplicates_rc(i,1);
                            multiple_match.guide_index_start(row,1) = num_duplicates_rc(i,1)+(sgRNA_len-1)+PAM_len;
                            multiple_match.dsb_index(row,1) = multiple_match.guide_index_start(row,1)-18-PAM_len;
                            multiple_match.dsb_coord(row,1) = reference_sequence.coord(multiple_match.dsb_index(row,1));
                            multiple_match.guide_coord_end(row,1) = reference_sequence.coord(num_duplicates_rc(i,1));
                            multiple_match.guide_coord_start(row,1) = reference_sequence.coord(num_duplicates_rc(i,1)+(sgRNA_len-1)+PAM_len);
                            row = row + 1;
                        end
                    elseif original_guides_full.category(mult_match_index(n,1)) == 2
                        num_duplicates = transpose(strfind(reference_seq,char(original_guides_full.guides_only(mult_match_index(n,1)))));
                        for i = 1:length(num_duplicates)
                            multiple_match.guides{row,1} = char(reference_seq((num_duplicates(i,1)-PAM_len):(num_duplicates(i,1)+(sgRNA_len-1))));
                            multiple_match.guides_only{row,1} = multiple_match.guides{row,1}((PAM_len+1):(PAM_len+sgRNA_len));
                            multiple_match.pam_only{row,1} = multiple_match.guides{row,1}(1:PAM_len);
                            multiple_match.pam(row,1) = {PAM_sequence1};
                            multiple_match.strand{row,1} = '+';
                            multiple_match.dsb_index(row,1) = num_duplicates(i,1)+18;
                            multiple_match.dsb_coord(row,1) = reference_sequence.coord(multiple_match.dsb_index(row,1));
                            multiple_match.guide_index_start(row,1) = num_duplicates(i,1)-PAM_len;
                            multiple_match.guide_index_end(row,1) = num_duplicates(i,1)+(sgRNA_len-1);
                            multiple_match.guide_coord_start(row,1) = reference_sequence.coord(num_duplicates(i,1)-PAM_len);
                            multiple_match.guide_coord_end(row,1) = reference_sequence.coord(num_duplicates(i,1)+(sgRNA_len-1));
                            row = row + 1;
                        end
                    elseif original_guides_full.category(mult_match_index(n,1)) == 3
                        num_duplicates_rc = transpose(strfind(reference_seq,seqrcomplement(char(original_guides_full.guides_only(mult_match_index(n,1))))));
                        for i = 1:length(num_duplicates_rc)
                            multiple_match.guides{row,1} = seqrcomplement(char(reference_seq(num_duplicates_rc(i,1):(num_duplicates_rc(i,1)+(sgRNA_len-1)+PAM_len))));
                            multiple_match.guides_only{row,1} = multiple_match.guides{row,1}(PAM_len+1:length(multiple_match.guides{row,1}));
                            multiple_match.pam_only{row,1} = multiple_match.guides{row,1}(1:PAM_len);
                            multiple_match.pam(row,1) = {PAM_sequence1};
                            multiple_match.strand{row,1} = '-';
                            multiple_match.guide_index_end(row,1) = num_duplicates_rc(i,1);
                            multiple_match.guide_index_start(row,1) = num_duplicates_rc(i,1)+(sgRNA_len-1)+PAM_len;
                            multiple_match.dsb_index(row,1) = multiple_match.guide_index_start(row,1)-18-PAM_len;
                            multiple_match.dsb_coord(row,1) = reference_sequence.coord(multiple_match.dsb_index(row,1));
                            multiple_match.guide_coord_end(row,1) = reference_sequence.coord(num_duplicates_rc(i,1));
                            multiple_match.guide_coord_start(row,1) = reference_sequence.coord(num_duplicates_rc(i,1)+(sgRNA_len-1)+PAM_len);
                            row = row + 1;
                        end
                    end
                    clearvars num_duplicates num_duplicates_rc
                end
            end
        end
        
        %% Haplotype Analysis
        
        if perform_haplotype_analysis == 1;
            disp('Haplotype Analysis...')
            temp_vcf_name  = [vcf_filenames(a,1),'.vcf'];
            identify_num_ind = 10000;
            if length(reference_seq)>algorithm_type
                fid = fopen(strjoin(temp_vcf_name,''));
                vcf_haplo_file5 = textscan(fid,'%s');
                fclose(fid);
                vcf_haplo_file4 = vcf_haplo_file5{1};
                clearvars vcf_haplo_file5
                string_identifier1 = '#CHROM';
                [answer_if1 vcf_haplo_file2] = ismember(string_identifier1,vcf_haplo_file4);
                vcf_haplo_file1 = vcf_haplo_file4(vcf_haplo_file2:length(vcf_haplo_file4));
                clearvars vcf_haplo_file4
                string_identifier2 = 'HG';
                string_identifier3 = 'NA';
                
                if length(vcf_haplo_file1)>identify_num_ind
                    vcf_haplo_file1_short = vcf_haplo_file1(1:identify_num_ind);
                else vcf_haplo_file1_short = vcf_haplo_file1;
                end
                
                num_individuals1 = length(cell2mat(strfind(vcf_haplo_file1_short,'HG')));
                num_individuals2 = length(cell2mat(strfind(vcf_haplo_file1_short,'NA')));
                num_individuals = num_individuals1+num_individuals2;
                num_haplo_vcf_columns = 9+num_individuals;
                clearvars vcf_haplo_file1_short
                if isempty(vcf_haplo_file1(num_haplo_vcf_columns+1:length(vcf_haplo_file1))) == 1
                    perform_haplotype_analysis = 0;
                end
                
                vcf_haplo_file = vcf_haplo_file1(num_haplo_vcf_columns+1:length(vcf_haplo_file1));
                indiv_ident.identifier = vcf_haplo_file1((num_haplo_vcf_columns-num_individuals+1):num_haplo_vcf_columns);
                indiv_ident.index = transpose(1:num_individuals);
                clearvars vcf_haplo_file1
                
                for h = 1:num_haplo_vcf_columns
                    dummy_var = h:num_haplo_vcf_columns:length(vcf_haplo_file);
                    eval(['column',num2str(h),'(:,1) = dummy_var;']);
                    clearvars dummy_var
                end
                vcf_haplo_snp.chr_id = str2double(vcf_haplo_file(column1,1));
                vcf_haplo_snp.bpstart = str2double(vcf_haplo_file(column2,1));
                vcf_haplo_snp.rs_id = vcf_haplo_file(column3,1);
                vcf_haplo_snp.ref_snp = vcf_haplo_file(column4,1);
                vcf_haplo_snp.alt_snp = vcf_haplo_file(column5,1);
                
                for i = 1:num_individuals;
                    eval(['vcf_haplo_snp.ind',num2str(i),'_full_genotype = vcf_haplo_file(column9+i,1);']);
                end
                total_length_haplo = length(vcf_haplo_file);
                num_rows_variants = (total_length_haplo/num_haplo_vcf_columns);
                clearvars vcf_haplo_file
                for i = 1:num_individuals
                    eval(['string_temp1 = char(vcf_haplo_snp.ind',num2str(i),'_full_genotype);']);
                    string_temp = string_temp1(:);
                    eval(['vcf_haplo_snp.ind',num2str(i),'_genotype1 = double(str2num(string_temp(1:num_rows_variants)));']);
                    eval(['vcf_haplo_snp.ind',num2str(i),'_genotype2 = double(str2num(string_temp(((2*num_rows_variants)+1):length(string_temp))));']);
                    clearvars string_temp string_temp1
                end
                
                for i = 1:num_haplo_vcf_columns
                    eval(['clearvars column',num2str(i),'']);
                end
                
            elseif length(reference_seq)<algorithm_type
                vcf_haplo_file4 = textread(char(strjoin(temp_vcf_name,'')),'%s');
                vcf_haplo_file3 = strfind(vcf_haplo_file4,'#CHROM');
                vcf_haplo_file2 = find(~cellfun(@isempty,vcf_haplo_file3));
                vcf_haplo_file1 = vcf_haplo_file4(vcf_haplo_file2:length(vcf_haplo_file4));
                num_individuals1 = length(cell2mat(strfind(vcf_haplo_file1,'HG')));
                num_individuals2 = length(cell2mat(strfind(vcf_haplo_file1,'NA')));
                num_individuals = num_individuals1+num_individuals2;
                num_haplo_vcf_columns = 9+num_individuals;
                if isempty(vcf_haplo_file1(num_haplo_vcf_columns+1:length(vcf_haplo_file1))) == 1
                    perform_haplotype_analysis = 0;
                end
                vcf_haplo_file = vcf_haplo_file1(num_haplo_vcf_columns+1:length(vcf_haplo_file1));
                indiv_ident.identifier = vcf_haplo_file1((num_haplo_vcf_columns-num_individuals+1):num_haplo_vcf_columns);
                indiv_ident.index = transpose(1:num_individuals);
                
                for n = 1:num_haplo_vcf_columns
                    eval(['row',num2str(n),'=',num2str(n),';']);
                    if n == num_haplo_vcf_columns
                        eval(['row',num2str(n+1),'=',num2str(n+1),';'])
                    end
                end
                
                for n = 1:length(vcf_haplo_file)
                    if row1 <= length(vcf_haplo_file)
                        eval(['vcf_haplo_snp.chr_id(n,1) = str2double(vcf_haplo_file(row1,1));']);
                        eval(['vcf_haplo_snp.bpstart(n,1) = str2double(vcf_haplo_file(row2,1));']);
                        eval(['vcf_haplo_snp.rs_id(n,1) = vcf_haplo_file(row3,1);']);
                        eval(['vcf_haplo_snp.ref_snp(n,1) = vcf_haplo_file(row4,1);']);
                        eval(['vcf_haplo_snp.alt_snp(n,1) = vcf_haplo_file(row5,1);']);
                        eval(['vcf_haplo_snp.qual(n,1) = vcf_haplo_file(row6,1);']);
                        eval(['vcf_haplo_snp.filter(n,1) = vcf_haplo_file(row7,1);']);
                        eval(['vcf_haplo_snp.info(n,1) = vcf_haplo_file(row8,1);']);
                        eval(['vcf_haplo_snp.format(n,1) = vcf_haplo_file(row9,1);']);
                        for i = 1:num_individuals;
                            eval(['vcf_haplo_snp.ind',num2str(i),'_full_genotype(n,1) = vcf_haplo_file(row9+i,1);']);
                            eval(['string_temp = char(vcf_haplo_snp.ind',num2str(i),'_full_genotype(n,1));']);
                            eval(['vcf_haplo_snp.ind',num2str(i),'_genotype1(n,1) = string_temp(1);']);
                            eval(['vcf_haplo_snp.ind',num2str(i),'_genotype2(n,1) = string_temp(3);']);
                        end
                        for j = 1:num_haplo_vcf_columns
                            eval(['row',num2str(j),' = row',num2str(j),'+num_haplo_vcf_columns;']);
                        end
                    elseif row1>length(vcf_haplo_file)
                        break
                    end
                end
                
                for i = 1:num_individuals
                    eval(['vcf_haplo_snp.ind',num2str(i),'_genotype1 = str2num(vcf_haplo_snp.ind',num2str(i),'_genotype1);']);
                    eval(['vcf_haplo_snp.ind',num2str(i),'_genotype2 = str2num(vcf_haplo_snp.ind',num2str(i),'_genotype2);']);
                end
            end
            
            num_haplotypes = 2*num_individuals;
            deletion = 0;
            
            for n = 1:num_haplotypes
                if mod(n,2) == 0
                    num_geno = 2;
                    num_ind = n-deletion;
                else
                    num_geno = 1;
                    num_ind = n-deletion;
                    deletion = deletion + 1;
                end
                cn_string = '<';
                eval(['haplotype_variants',num2str(n),'.index = find(vcf_haplo_snp.ind',num2str(num_ind),'_genotype',num2str(num_geno),' > 0);']);
                eval(['haplotype_variants',num2str(n),'.chr_id = vcf_haplo_snp.chr_id(haplotype_variants',num2str(n),'.index);'])
                eval(['haplotype_variants',num2str(n),'.variant_coord_start = vcf_haplo_snp.bpstart(haplotype_variants',num2str(n),'.index);'])
                eval(['haplotype_variants',num2str(n),'.rs_id = vcf_haplo_snp.rs_id(haplotype_variants',num2str(n),'.index);'])
                eval(['haplotype_variants',num2str(n),'.ref_snp = vcf_haplo_snp.ref_snp(haplotype_variants',num2str(n),'.index);'])
                eval(['haplotype_variants',num2str(n),'.alt_snp = vcf_haplo_snp.alt_snp(haplotype_variants',num2str(n),'.index);'])
                eval(['haplotype_variants',num2str(n),'.genotype = vcf_haplo_snp.ind',num2str(num_ind),'_genotype',num2str(num_geno),'(haplotype_variants',num2str(n),'.index);'])
                eval(['[variant_index_start_temp haplotype_variants',num2str(n),'.variant_index_start] = ismember(haplotype_variants',num2str(n),'.variant_coord_start,reference_sequence.coord);'])
                eval(['test_for_cn = vcf_haplo_snp.alt_snp(haplotype_variants',num2str(n),'.index);'])
                if sum(sum(ismember(char(test_for_cn),cn_string),2))>0
                    eval(['haplotype_variants',num2str(n),'.index = find(sum(ismember(char(haplotype_variants',num2str(n),'.alt_snp),cn_string),2)<1);']);
                    eval(['haplotype_variants',num2str(n),'.chr_id = haplotype_variants',num2str(n),'.chr_id(haplotype_variants',num2str(n),'.index);'])
                    eval(['haplotype_variants',num2str(n),'.variant_coord_start = haplotype_variants',num2str(n),'.variant_coord_start(haplotype_variants',num2str(n),'.index);'])
                    eval(['haplotype_variants',num2str(n),'.rs_id = haplotype_variants',num2str(n),'.rs_id(haplotype_variants',num2str(n),'.index);'])
                    eval(['haplotype_variants',num2str(n),'.ref_snp = haplotype_variants',num2str(n),'.ref_snp(haplotype_variants',num2str(n),'.index);'])
                    eval(['haplotype_variants',num2str(n),'.alt_snp = haplotype_variants',num2str(n),'.alt_snp(haplotype_variants',num2str(n),'.index);'])
                    eval(['haplotype_variants',num2str(n),'.genotype = haplotype_variants',num2str(n),'.genotype(haplotype_variants',num2str(n),'.index);'])
                    eval(['haplotype_variants',num2str(n),'.variant_index_start = haplotype_variants',num2str(n),'.variant_index_start(haplotype_variants',num2str(n),'.index);'])
                end
            end
            
            row_sub = 1;
            row_del = 1;
            row_ins = 1;
            disp('Classifying Variants...')
            if length(reference_seq)<algorithm_type
                for n = 1:num_haplotypes
                    nucleotide_list = ['AGCT'];
                    eval(['num_var_haplo_iterate = length(haplotype_variants',num2str(n),'.ref_snp);'])
                    for i = 1:num_var_haplo_iterate
                        eval(['temp1 = length(char(haplotype_variants',num2str(n),'.ref_snp(i,1)));']);
                        eval(['temp2 = length(char(haplotype_variants',num2str(n),'.alt_snp(i,1)));']);
                        eval(['temp3 = char(haplotype_variants',num2str(n),'.alt_snp(i,1));']);
                        
                        alt_snp_complex = 0;
                        if isempty(strfind(temp3,',')) == 0;
                            if isempty(strfind(temp3,'"')) == 0
                                temp3 = temp3(2:(length(temp3)-1));
                            end
                            temp_list = transpose(strsplit(temp3,','));
                            eval(['haplotype_variants',num2str(n),'.alt_snp(i,1) = temp_list(haplotype_variants',num2str(n),'.genotype(i,1));']);
                            eval(['temp2 = length(char(haplotype_variants',num2str(n),'.alt_snp(i,1)));']);
                            alt_snp_complex = 1;
                            
                        end
                        eval(['haplotype_variants',num2str(n),'.complex(i,1) = alt_snp_complex;']);
                        for g = 1:length(temp2)
                            if temp1(g,1)==temp2(g,1)
                                eval(['haplotype_variants',num2str(n),'.identifier(i,1) = 1;']);
                                eval(['haplotype_variants',num2str(n),'.variant_index_start(i,1) = find(haplotype_variants',num2str(n),'.variant_coord_start(i,1) == reference_sequence.coord);'])
                                eval(['haplotype_variants',num2str(n),'.variant_index_end(i,1) = haplotype_variants',num2str(n),'.variant_index_start(i,1);']);
                                eval(['haplotype_variants',num2str(n),'.variant_coord_end(i,1) = haplotype_variants',num2str(n),'.variant_coord_start(i,1);']);
                                eval(['haplotype_variants',num2str(n),'.variant_len(i,1) = 0;']);
                                eval(['substitutions.bpstart(row_sub,1) = haplotype_variants',num2str(n),'.variant_coord_start(i,1);']);
                                eval(['substitutions.bpend(row_sub,1) = haplotype_variants',num2str(n),'.variant_coord_end(i,1);']);
                                eval(['substitutions.identifier(row_sub,1) = haplotype_variants',num2str(n),'.identifier(i,1);']);
                                row_sub = row_sub + 1;
                                if g > 1
                                    haplo_sub = '1';
                                    eval(['haplotype_variants',num2str(n),'.identifier(i,1) = strcat(num2str(haplotype_variants',num2str(n),'.identifier(i,1),haplo_sub);']);
                                end
                            elseif temp1(g,1)>temp2(g,1)
                                eval(['haplotype_variants',num2str(n),'.identifier(i,1) = 2;']);
                                eval(['haplotype_variants',num2str(n),'.variant_len(i,1) = size(char(haplotype_variants',num2str(n),'.ref_snp(i,1)),2)-1;'])
                                eval(['haplotype_variants',num2str(n),'.variant_index_start(i,1) = find(haplotype_variants',num2str(n),'.variant_coord_start(i,1) == reference_sequence.coord);'])
                                eval(['haplotype_variants',num2str(n),'.variant_index_end(i,1) = haplotype_variants',num2str(n),'.variant_index_start(i,1)+haplotype_variants',num2str(n),'.variant_len(i,1);'])
                                eval(['haplotype_variants',num2str(n),'.variant_coord_end(i,1) = haplotype_variants',num2str(n),'.variant_coord_start(i,1)+haplotype_variants',num2str(n),'.variant_len(i,1);'])
                                eval(['deletions.bpstart(row_del,1) = haplotype_variants',num2str(n),'.variant_coord_start(i,1);']);
                                eval(['deletions.bpend(row_del,1) = haplotype_variants',num2str(n),'.variant_coord_end(i,1);']);
                                eval(['deletions.identifier(row_del,1) = haplotype_variants',num2str(n),'.identifier(i,1);']);
                                row_del = row_del + 1;
                                if g > 1
                                    haplo_del = '2';
                                    eval(['haplotype_variants',num2str(n),'.identifier(i,1) = strcat(num2str(haplotype_variants',num2str(n),'.identifier(i,1),haplo_del);']);
                                end
                            elseif temp1(g,1)<temp2(g,1)
                                eval(['haplotype_variants',num2str(n),'.identifier(i,1) = 3;']);
                                eval(['haplotype_variants',num2str(n),'.variant_len(i,1) = size(char(haplotype_variants',num2str(n),'.alt_snp(i,1)),2)-1;'])
                                eval(['haplotype_variants',num2str(n),'.variant_index_start(i,1) = find(haplotype_variants',num2str(n),'.variant_coord_start(i,1) == reference_sequence.coord);'])
                                eval(['haplotype_variants',num2str(n),'.variant_index_end(i,1) = haplotype_variants',num2str(n),'.variant_index_start(i,1)+haplotype_variants',num2str(n),'.variant_len(i,1);'])
                                eval(['haplotype_variants',num2str(n),'.variant_coord_end(i,1) = haplotype_variants',num2str(n),'.variant_coord_start(i,1)+haplotype_variants',num2str(n),'.variant_len(i,1);'])
                                eval(['insertions.bpstart(row_ins,1) = haplotype_variants',num2str(n),'.variant_coord_start(i,1);']);
                                eval(['insertions.bpend(row_ins,1) = haplotype_variants',num2str(n),'.variant_coord_end(i,1);']);
                                eval(['insertions.identifier(row_ins,1) = haplotype_variants',num2str(n),'.identifier(i,1);']);
                                row_ins = row_ins + 1;
                                if g > 1
                                    haplo_ins = '3';
                                    eval(['haplotype_variants',num2str(n),'.identifier(i,1) = strcat(num2str(haplotype_variants',num2str(n),'.identifier(i,1),haplo_ins);']);
                                end
                                
                            end
                        end
                    end
                end
                
            elseif length(reference_seq)>algorithm_type
                for n = 1:num_haplotypes
                    eval(['temp1 = cellfun(''length'',haplotype_variants',num2str(n),'.ref_snp);']);
                    alt_snp_complex = 0;
                    if eval(['ismember('','',char(haplotype_variants',num2str(n),'.alt_snp)) == 1']);
                        eval(['temp3 = char(haplotype_variants',num2str(n),'.alt_snp);']);
                        eval(['temp4 = haplotype_variants',num2str(n),'.alt_snp;']);
                        eval(['temp3_index2 = strfind(haplotype_variants',num2str(n),'.alt_snp,'','');'])
                        temp3_index1 = cellfun('length',temp3_index2);
                        temp3_index = find(temp3_index1>0);
                        for y = 1:length(temp3_index);
                            if ismember('"',temp3(temp3_index(y,1))) == 1
                                temp5 = char(temp4(temp3_index(y,1)));
                                temp4{temp3_index(y,1)} = temp5(2:(length(temp5)-1));
                            end
                            temp_list = transpose(strsplit(char(temp4(temp3_index(y,1))),','));
                            eval(['haplotype_variants',num2str(n),'.alt_snp(temp3_index(y,1)) = temp_list(haplotype_variants',num2str(n),'.genotype(temp3_index(y,1)));']);
                            alt_snp_complex = 1;
                        end
                        
                    end
                    
                    eval(['temp2 = cellfun(''length'',haplotype_variants',num2str(n),'.alt_snp);']);
                    if exist('temp1','var') == 1
                        temp_diff = temp2-temp1;
                        temp_sub_index = find(temp_diff == 0);
                        temp_del_index = find(temp_diff < 0);
                        temp_ins_index = find(temp_diff > 0);
                        
                        if length(temp_sub_index) > 0
                            eval(['haplotype_variants',num2str(n),'.identifier(temp_sub_index,1) = 1;']);
                            eval(['haplotype_variants',num2str(n),'.variant_index_end(temp_sub_index,1) = haplotype_variants',num2str(n),'.variant_index_start(temp_sub_index,1);']);
                            eval(['haplotype_variants',num2str(n),'.variant_coord_end(temp_sub_index,1) = haplotype_variants',num2str(n),'.variant_coord_start(temp_sub_index,1);']);
                            eval(['haplotype_variants',num2str(n),'.variant_len(temp_sub_index,1) = 0;']);
                            eval(['substitutions.bpstart = haplotype_variants',num2str(n),'.variant_coord_start(temp_sub_index,1);']);
                            eval(['substitutions.bpend = haplotype_variants',num2str(n),'.variant_coord_end(temp_sub_index,1);']);
                            eval(['substitutions.identifier = haplotype_variants',num2str(n),'.identifier(temp_sub_index,1);']);
                        end
                        
                        if length(temp_del_index) > 0
                            eval(['haplotype_variants',num2str(n),'.identifier(temp_del_index,1) = 2;']);
                            eval(['haplotype_variants',num2str(n),'.variant_len(temp_del_index,1) = abs(temp_diff(temp_del_index))-1;'])
                            eval(['haplotype_variants',num2str(n),'.variant_index_end(temp_del_index,1) = haplotype_variants',num2str(n),'.variant_index_start(temp_del_index,1)+(abs(temp_diff(temp_del_index)));'])
                            eval(['haplotype_variants',num2str(n),'.variant_coord_end(temp_del_index,1) = haplotype_variants',num2str(n),'.variant_coord_start(temp_del_index,1)+(abs(temp_diff(temp_del_index)));'])
                            eval(['deletions.bpstart = haplotype_variants',num2str(n),'.variant_coord_start(temp_del_index,1);']);
                            eval(['deletions.bpend = haplotype_variants',num2str(n),'.variant_coord_end(temp_del_index,1);']);
                            eval(['deletions.identifier = haplotype_variants',num2str(n),'.identifier(temp_del_index,1);']);
                        end
                        
                        if length(temp_ins_index) > 0
                            eval(['haplotype_variants',num2str(n),'.identifier(temp_ins_index,1) = 3;']);
                            eval(['haplotype_variants',num2str(n),'.variant_len(temp_ins_index,1) = temp_diff(temp_ins_index)-1;'])
                            eval(['haplotype_variants',num2str(n),'.variant_index_end(temp_ins_index,1) = haplotype_variants',num2str(n),'.variant_index_start(temp_ins_index,1)+(abs(temp_diff(temp_ins_index))-1);'])
                            eval(['haplotype_variants',num2str(n),'.variant_coord_end(temp_ins_index,1) = haplotype_variants',num2str(n),'.variant_coord_start(temp_ins_index,1)+(abs(temp_diff(temp_ins_index))-1);'])
                            eval(['insertions.bpstart = haplotype_variants',num2str(n),'.variant_coord_start(temp_ins_index,1);']);
                            eval(['insertions.bpend = haplotype_variants',num2str(n),'.variant_coord_end(temp_ins_index,1);']);
                            eval(['insertions.identifier = haplotype_variants',num2str(n),'.identifier(temp_ins_index,1);']);
                        end
                        clearvars temp1 temp2 temp3 temp4 temp5 temp3_index temp3_index1 temp3_index2 temp_sub_index temp_del_index temp_ins_index
                    end
                end
            end
            
            determine_hap_number = num_haplotypes;
            for n = 1:determine_hap_number
                if eval(['isempty(haplotype_variants',num2str(n),'.index) == 0'])
                    if exist('num_haplotypes_list','var') == 0;
                        num_haplotypes_list = n;
                    else
                        num_haplotypes_list = [num_haplotypes_list;n];
                    end
                end
            end
            
            for n = 1:num_haplotypes
                if isempty(find(num_haplotypes_list==n)) == 0
                    eval(['[haplotype_variants',num2str(n),'.identifier haplo_index] = sort(haplotype_variants',num2str(n),'.identifier,1);'])
                    eval(['haplotype_variants',num2str(n),'.chr_id = haplotype_variants',num2str(n),'.chr_id(haplo_index);']);
                    eval(['haplotype_variants',num2str(n),'.variant_index_start = haplotype_variants',num2str(n),'.variant_index_start(haplo_index);']);
                    eval(['haplotype_variants',num2str(n),'.variant_index_end = haplotype_variants',num2str(n),'.variant_index_end(haplo_index);']);
                    eval(['haplotype_variants',num2str(n),'.rs_id = haplotype_variants',num2str(n),'.rs_id(haplo_index);']);
                    eval(['haplotype_variants',num2str(n),'.variant_coord_start = haplotype_variants',num2str(n),'.variant_coord_start(haplo_index);']);
                    eval(['haplotype_variants',num2str(n),'.variant_coord_end = haplotype_variants',num2str(n),'.variant_coord_end(haplo_index);']);
                    eval(['haplotype_variants',num2str(n),'.ref_snp = haplotype_variants',num2str(n),'.ref_snp(haplo_index);']);
                    eval(['haplotype_variants',num2str(n),'.alt_snp = haplotype_variants',num2str(n),'.alt_snp(haplo_index);']);
                    eval(['haplotype_variants',num2str(n),'.variant_len = haplotype_variants',num2str(n),'.variant_len(haplo_index);']);
                    eval(['haplotype_variants',num2str(n),'.genotype = haplotype_variants',num2str(n),'.genotype(haplo_index);']);
                end
            end
            
            full_list_haplo_variants.chr_id = vcf_haplo_snp.chr_id;
            full_list_haplo_variants.ref_snp = vcf_haplo_snp.ref_snp;
            full_list_haplo_variants.alt_snp = vcf_haplo_snp.alt_snp;
            full_list_haplo_variants.rs_id = vcf_haplo_snp.rs_id;
            full_list_haplo_variants.bpstart = vcf_haplo_snp.bpstart;
            for g = 1:length(full_list_haplo_variants.bpstart)
                full_list_haplo_variants.variant_index_start(g,1) = find(full_list_haplo_variants.bpstart(g,1)==reference_sequence.coord);
            end
            
            eval(['temp1 = cellfun(''length'',full_list_haplo_variants.ref_snp);']);
            eval(['temp2 = cellfun(''length'',full_list_haplo_variants.alt_snp);']);
            temp_diff = temp2-temp1;
            temp_sub_index = find(temp_diff == 0);
            temp_del_index = find(temp_diff < 0);
            temp_ins_index = find(temp_diff > 0);
            full_list_haplo_variants.identifier(1:length(full_list_haplo_variants.bpstart),1) = {'Identifier'};
            full_list_haplo_variants.identifier(temp_sub_index,1) = {'Substitution'};
            full_list_haplo_variants.identifier(temp_del_index,1) = {'Deletion'};
            full_list_haplo_variants.identifier(temp_ins_index,1) = {'Insertion'};
            full_list_haplo_variants.variant_len(1:length(full_list_haplo_variants.bpstart),1) = 0;
            full_list_haplo_variants.variant_len(temp_sub_index,1) = 0;
            full_list_haplo_variants.variant_len(temp_del_index,1) = abs(temp_diff(temp_del_index))-1;
            full_list_haplo_variants.variant_len(temp_ins_index,1) = temp_diff(temp_ins_index)-1;
            full_list_haplo_variants.bpend = full_list_haplo_variants.bpstart+full_list_haplo_variants.variant_len;
            full_list_haplo_variants.variant_index_end = full_list_haplo_variants.variant_index_start+full_list_haplo_variants.variant_len;
            
            disp('Building Haplotype Sequences...')
            
            if length(reference_seq) < algorithm_type2
                for n = 1:num_haplotypes
                    if isempty(find(num_haplotypes_list==n)) == 0
                        row_w = 1;
                        eval(['num_variants_iterate = length(haplotype_variants',num2str(n),'.ref_snp);']);
                        full_sequence_modify.nucleotide = reference_sequence.nucleotide;
                        full_sequence_modify.index = reference_sequence.number;
                        full_sequence_modify.coord = reference_sequence.coord;
                        for j = 1:num_variants_iterate;
                            if eval(['haplotype_variants',num2str(n),'.identifier(j,1) == 1']) & eval(['strcmp(full_sequence_modify.nucleotide(find(haplotype_variants',num2str(n),'.variant_index_start(j,1)==full_sequence_modify.index)),haplotype_variants',num2str(n),'.ref_snp(j,1)) == 1']);
                                eval(['haplotype',num2str(n),'.ref_snp(row_w,1) = haplotype_variants',num2str(n),'.ref_snp(j,1);']);
                                eval(['haplotype',num2str(n),'.alt_snp(row_w,1) = haplotype_variants',num2str(n),'.alt_snp(j,1);']);
                                eval(['haplotype',num2str(n),'.from_seq(row_w,1) = full_sequence_modify.nucleotide(find(haplotype_variants',num2str(n),'.variant_index_start(j,1)==full_sequence_modify.index));']);
                                eval(['haplotype',num2str(n),'.variant_index_start(row_w,1) = haplotype_variants',num2str(n),'.variant_index_start(j,1);']);
                                eval(['haplotype',num2str(n),'.variant_index_end(row_w,1) = haplotype_variants',num2str(n),'.variant_index_end(j,1);']);
                                eval(['haplotype',num2str(n),'.variant_coord_start(row_w,1) = haplotype_variants',num2str(n),'.variant_coord_start(j,1);']);
                                eval(['haplotype',num2str(n),'.variant_coord_end(row_w,1) = haplotype_variants',num2str(n),'.variant_coord_end(j,1);']);
                                eval(['haplotype',num2str(n),'.rs_id(row_w,1) = haplotype_variants',num2str(n),'.rs_id(j,1);']);
                                eval(['haplotype',num2str(n),'.identifier(row_w,1) = haplotype_variants',num2str(n),'.identifier(j,1);']);
                                eval(['haplotype',num2str(n),'.variant_len(row_w,1) = haplotype_variants',num2str(n),'.variant_len(j,1);']);
                                eval(['haplotype',num2str(n),'.type{row_w,1} = ''Substitution'';']);
                                eval(['full_sequence_modify.nucleotide(haplotype_variants',num2str(n),'.variant_index_start(j,1)) = haplotype_variants',num2str(n),'.alt_snp(j,1);'])
                                row_w = row_w + 1;
                            elseif eval(['haplotype_variants',num2str(n),'.identifier(j,1) == 2']) & eval(['strcmp(strjoin(full_sequence_modify.nucleotide(haplotype_variants',num2str(n),'.variant_index_start(j,1):haplotype_variants',num2str(n),'.variant_index_end(j,1)),''''),haplotype_variants',num2str(n),'.ref_snp(j,1)) == 1']);
                                eval(['haplotype',num2str(n),'.ref_snp(row_w,1) = haplotype_variants',num2str(n),'.ref_snp(j,1);']);
                                eval(['haplotype',num2str(n),'.alt_snp(row_w,1) = haplotype_variants',num2str(n),'.alt_snp(j,1);']);
                                eval(['haplotype',num2str(n),'.from_seq(row_w,1) = full_sequence_modify.nucleotide(find(haplotype_variants',num2str(n),'.variant_index_start(j,1)==full_sequence_modify.index));']);
                                eval(['haplotype',num2str(n),'.variant_index_start(row_w,1) = haplotype_variants',num2str(n),'.variant_index_start(j,1);']);
                                eval(['haplotype',num2str(n),'.variant_index_end(row_w,1) = haplotype_variants',num2str(n),'.variant_index_end(j,1);']);
                                eval(['haplotype',num2str(n),'.variant_coord_start(row_w,1) = haplotype_variants',num2str(n),'.variant_coord_start(j,1);']);
                                eval(['haplotype',num2str(n),'.variant_coord_end(row_w,1) = haplotype_variants',num2str(n),'.variant_coord_end(j,1);']);
                                eval(['haplotype',num2str(n),'.rs_id(row_w,1) = haplotype_variants',num2str(n),'.rs_id(j,1);']);
                                eval(['haplotype',num2str(n),'.identifier(row_w,1) = haplotype_variants',num2str(n),'.identifier(j,1);']);
                                eval(['haplotype',num2str(n),'.variant_len(row_w,1) = haplotype_variants',num2str(n),'.variant_len(j,1);']);
                                eval(['haplotype',num2str(n),'.type{row_w,1} = ''Deletion'';']);
                                eval(['full_sequence_modify.nucleotide((haplotype_variants',num2str(n),'.variant_index_start(j,1)+(length(char(haplotype_variants',num2str(n),'.alt_snp(j,1))))):haplotype_variants',num2str(n),'.variant_index_end(j,1)) = {''D''};']);
                                row_w = row_w + 1;
                            elseif eval(['haplotype_variants',num2str(n),'.identifier(j,1) == 3'])
                                eval(['haplotype',num2str(n),'.ref_snp(row_w,1) = haplotype_variants',num2str(n),'.ref_snp(j,1);']);
                                eval(['haplotype',num2str(n),'.alt_snp(row_w,1) = haplotype_variants',num2str(n),'.alt_snp(j,1);']);
                                eval(['haplotype',num2str(n),'.variant_index_start(row_w,1) = haplotype_variants',num2str(n),'.variant_index_start(j,1);']);
                                eval(['haplotype',num2str(n),'.variant_index_end(row_w,1) = haplotype_variants',num2str(n),'.variant_index_end(j,1);']);
                                eval(['haplotype',num2str(n),'.variant_coord_start(row_w,1) = haplotype_variants',num2str(n),'.variant_coord_start(j,1);']);
                                eval(['haplotype',num2str(n),'.variant_coord_end(row_w,1) = haplotype_variants',num2str(n),'.variant_coord_end(j,1);']);
                                eval(['haplotype',num2str(n),'.from_seq(row_w,1) = full_sequence_modify.nucleotide(find(haplotype_variants',num2str(n),'.variant_index_start(j,1)==full_sequence_modify.index));']);
                                eval(['haplotype',num2str(n),'.rs_id(row_w,1) = haplotype_variants',num2str(n),'.rs_id(j,1);']);
                                eval(['haplotype',num2str(n),'.identifier(row_w,1) = haplotype_variants',num2str(n),'.identifier(j,1);']);
                                eval(['haplotype',num2str(n),'.variant_len(row_w,1) = haplotype_variants',num2str(n),'.variant_len(j,1);']);
                                eval(['haplotype',num2str(n),'.type{row_w,1} = ''Insertion'';']);
                                ins_length_increase = 0;
                                eval(['temp_alt_snp = haplotype_variants',num2str(n),'.alt_snp(j,1);;'])
                                split_ins = transpose(cellstr(char(temp_alt_snp)')');
                                split_ins = split_ins(2:length(split_ins));
                                ins_length = length(split_ins);
                                eval(['index_ins_original = find(haplotype_variants',num2str(n),'.variant_index_start(j,1)==full_sequence_modify.index);']);
                                eval(['index_ins = find(haplotype_variants',num2str(n),'.variant_index_start(j,1)==full_sequence_modify.index);']);
                                eval(['ref_len = length(char(haplotype_variants',num2str(n),'.ref_snp(j,1)))-1;'])
                                eval(['alt_len = length(char(haplotype_variants',num2str(n),'.alt_snp(j,1)))-1;'])
                                diff_factor = alt_len-ref_len;
                                full_sequence_modify1 = full_sequence_modify.nucleotide(1:(index_ins+ref_len));
                                full_sequence_modify1((index_ins+ref_len+diff_factor+1):(length(full_sequence_modify.nucleotide)+diff_factor)) = full_sequence_modify.nucleotide(index_ins+ref_len+1:length(full_sequence_modify.nucleotide));
                                full_sequence_modify.nucleotide = full_sequence_modify1;
                                full_sequence_modify_index1 = full_sequence_modify.index(1:(index_ins+ref_len));
                                full_sequence_modify_index1((index_ins+ref_len+diff_factor+1):(length(full_sequence_modify.index)+diff_factor)) = full_sequence_modify.index(index_ins+ref_len+1:length(full_sequence_modify.index));
                                full_sequence_modify.index = full_sequence_modify_index1;
                                full_sequence_modify_coord1 = full_sequence_modify.coord(1:(index_ins+ref_len));
                                full_sequence_modify_coord1((index_ins+ref_len+diff_factor+1):(length(full_sequence_modify.coord)+diff_factor)) = full_sequence_modify.coord(index_ins+ref_len+1:length(full_sequence_modify.coord));
                                full_sequence_modify.coord = full_sequence_modify_coord1;
                                
                                index_replace = index_ins+1;
                                ins_count = 0;
                                for i = 1:ins_length
                                    full_sequence_modify.nucleotide(index_replace+ins_count) = split_ins(i,1);
                                    ins_count = ins_count + 1;
                                end
                                
                                row_w = row_w + 1;
                            end
                        end
                        
                        ident_del_sites = strfind(full_sequence_modify.nucleotide,'D');
                        ident_del_sites = find(~cellfun(@isempty,ident_del_sites));
                        full_sequence_modify.nucleotide(ident_del_sites) = [];
                        eval(['haplotype',num2str(n),'.sequence = strjoin(transpose(full_sequence_modify.nucleotide),'''');']);
                    end
                end
                
                clearvars full_sequence_modify
            elseif length(reference_seq) > algorithm_type2
                for n = 1:num_haplotypes
                    if isempty(find(num_haplotypes_list==n)) == 0
                        full_sequence_modify.nucleotide = reference_sequence.nucleotide;
                        full_sequence_modify.index = reference_sequence.number;
                        full_sequence_modify.coord = reference_sequence.coord;
                        identifier_sub_string = 'Substitution';
                        if eval(['length(find(haplotype_variants',num2str(n),'.identifier == 1))>0']);
                            eval(['haplo_sub_index = find(haplotype_variants',num2str(n),'.identifier == 1);']);;
                            eval(['haplotype',num2str(n),'.ref_snp(haplo_sub_index,1) = haplotype_variants',num2str(n),'.ref_snp(haplo_sub_index,1);']);
                            eval(['haplotype',num2str(n),'.alt_snp(haplo_sub_index,1) = haplotype_variants',num2str(n),'.alt_snp(haplo_sub_index,1);']);
                            eval(['haplotype',num2str(n),'.from_seq(haplo_sub_index,1) = full_sequence_modify.nucleotide(haplotype_variants',num2str(n),'.variant_index_start(haplo_sub_index,1));']);
                            eval(['haplotype',num2str(n),'.variant_index_start(haplo_sub_index,1) = haplotype_variants',num2str(n),'.variant_index_start(haplo_sub_index,1);']);
                            eval(['haplotype',num2str(n),'.variant_index_end(haplo_sub_index,1) = haplotype_variants',num2str(n),'.variant_index_end(haplo_sub_index,1);']);
                            eval(['haplotype',num2str(n),'.variant_coord_start(haplo_sub_index,1) = haplotype_variants',num2str(n),'.variant_coord_start(haplo_sub_index,1);']);
                            eval(['haplotype',num2str(n),'.variant_coord_end(haplo_sub_index,1) = haplotype_variants',num2str(n),'.variant_coord_end(haplo_sub_index,1);']);
                            eval(['haplotype',num2str(n),'.rs_id(haplo_sub_index,1) = haplotype_variants',num2str(n),'.rs_id(haplo_sub_index,1);']);
                            eval(['haplotype',num2str(n),'.identifier(haplo_sub_index,1) = haplotype_variants',num2str(n),'.identifier(haplo_sub_index,1);']);
                            eval(['haplotype',num2str(n),'.variant_len(haplo_sub_index,1) = haplotype_variants',num2str(n),'.variant_len(haplo_sub_index,1);']);
                            eval(['full_sequence_modify.nucleotide(haplotype_variants',num2str(n),'.variant_index_start(haplo_sub_index,1)) = haplotype_variants',num2str(n),'.alt_snp(haplo_sub_index,1);'])
                        end
                        eval(['num_variants_iterate = length(haplotype_variants',num2str(n),'.ref_snp);']);
                        if exist('haplo_sub_index','var') == 1
                            row_w = length(haplo_sub_index) + 1;
                        elseif exist('haplo_sub_index','var') == 0
                            row_w = 1;
                        end
                        for j = 1:num_variants_iterate;
                            if eval(['haplotype_variants',num2str(n),'.identifier(j,1) == 2']) & eval(['strcmp(strjoin(full_sequence_modify.nucleotide(haplotype_variants',num2str(n),'.variant_index_start(j,1):haplotype_variants',num2str(n),'.variant_index_end(j,1)),''''),haplotype_variants',num2str(n),'.ref_snp(j,1)) == 1']);
                                eval(['haplotype',num2str(n),'.ref_snp(row_w,1) = haplotype_variants',num2str(n),'.ref_snp(j,1);']);
                                eval(['haplotype',num2str(n),'.alt_snp(row_w,1) = haplotype_variants',num2str(n),'.alt_snp(j,1);']);
                                eval(['haplotype',num2str(n),'.from_seq(row_w,1) = full_sequence_modify.nucleotide(find(haplotype_variants',num2str(n),'.variant_index_start(j,1)==full_sequence_modify.index));']);
                                eval(['haplotype',num2str(n),'.variant_index_start(row_w,1) = haplotype_variants',num2str(n),'.variant_index_start(j,1);']);
                                eval(['haplotype',num2str(n),'.variant_index_end(row_w,1) = haplotype_variants',num2str(n),'.variant_index_end(j,1);']);
                                eval(['haplotype',num2str(n),'.variant_coord_start(row_w,1) = haplotype_variants',num2str(n),'.variant_coord_start(j,1);']);
                                eval(['haplotype',num2str(n),'.variant_coord_end(row_w,1) = haplotype_variants',num2str(n),'.variant_coord_end(j,1);']);
                                eval(['haplotype',num2str(n),'.rs_id(row_w,1) = haplotype_variants',num2str(n),'.rs_id(j,1);']);
                                eval(['haplotype',num2str(n),'.identifier(row_w,1) = haplotype_variants',num2str(n),'.identifier(j,1);']);
                                eval(['haplotype',num2str(n),'.variant_len(row_w,1) = haplotype_variants',num2str(n),'.variant_len(j,1);']);
                                eval(['full_sequence_modify.nucleotide((haplotype_variants',num2str(n),'.variant_index_start(j,1)+(length(char(haplotype_variants',num2str(n),'.alt_snp(j,1))))):haplotype_variants',num2str(n),'.variant_index_end(j,1)) = {''D''};']);
                                row_w = row_w + 1;
                            elseif eval(['haplotype_variants',num2str(n),'.identifier(j,1) == 3'])
                                eval(['haplotype',num2str(n),'.ref_snp(row_w,1) = haplotype_variants',num2str(n),'.ref_snp(j,1);']);
                                eval(['haplotype',num2str(n),'.alt_snp(row_w,1) = haplotype_variants',num2str(n),'.alt_snp(j,1);']);
                                eval(['haplotype',num2str(n),'.variant_index_start(row_w,1) = haplotype_variants',num2str(n),'.variant_index_start(j,1);']);
                                eval(['haplotype',num2str(n),'.variant_index_end(row_w,1) = haplotype_variants',num2str(n),'.variant_index_end(j,1);']);
                                eval(['haplotype',num2str(n),'.variant_coord_start(row_w,1) = haplotype_variants',num2str(n),'.variant_coord_start(j,1);']);
                                eval(['haplotype',num2str(n),'.variant_coord_end(row_w,1) = haplotype_variants',num2str(n),'.variant_coord_end(j,1);']);
                                eval(['haplotype',num2str(n),'.from_seq(row_w,1) = full_sequence_modify.nucleotide(find(haplotype_variants',num2str(n),'.variant_index_start(j,1)==full_sequence_modify.index));']);
                                eval(['haplotype',num2str(n),'.rs_id(row_w,1) = haplotype_variants',num2str(n),'.rs_id(j,1);']);
                                eval(['haplotype',num2str(n),'.identifier(row_w,1) = haplotype_variants',num2str(n),'.identifier(j,1);']);
                                eval(['haplotype',num2str(n),'.variant_len(row_w,1) = haplotype_variants',num2str(n),'.variant_len(j,1);']);
                                ins_length_increase = 0;
                                eval(['temp_alt_snp = haplotype_variants',num2str(n),'.alt_snp(j,1);;'])
                                split_ins = transpose(cellstr(char(temp_alt_snp)')');
                                split_ins = split_ins(2:length(split_ins));
                                ins_length = length(split_ins);
                                eval(['index_ins_original = find(haplotype_variants',num2str(n),'.variant_index_start(j,1)==full_sequence_modify.index);']);
                                eval(['index_ins = find(haplotype_variants',num2str(n),'.variant_index_start(j,1)==full_sequence_modify.index);']);
                                eval(['ref_len = length(char(haplotype_variants',num2str(n),'.ref_snp(j,1)))-1;'])
                                eval(['alt_len = length(char(haplotype_variants',num2str(n),'.alt_snp(j,1)))-1;'])
                                diff_factor = alt_len-ref_len;
                                full_sequence_modify1 = full_sequence_modify.nucleotide(1:(index_ins+ref_len));
                                full_sequence_modify1((index_ins+ref_len+diff_factor+1):(length(full_sequence_modify.nucleotide)+diff_factor)) = full_sequence_modify.nucleotide(index_ins+ref_len+1:length(full_sequence_modify.nucleotide));
                                full_sequence_modify.nucleotide = full_sequence_modify1;
                                full_sequence_modify_index1 = full_sequence_modify.index(1:(index_ins+ref_len));
                                full_sequence_modify_index1((index_ins+ref_len+diff_factor+1):(length(full_sequence_modify.index)+diff_factor)) = full_sequence_modify.index(index_ins+ref_len+1:length(full_sequence_modify.index));
                                full_sequence_modify.index = full_sequence_modify_index1;
                                full_sequence_modify_coord1 = full_sequence_modify.coord(1:(index_ins+ref_len));
                                full_sequence_modify_coord1((index_ins+ref_len+diff_factor+1):(length(full_sequence_modify.coord)+diff_factor)) = full_sequence_modify.coord(index_ins+ref_len+1:length(full_sequence_modify.coord));
                                full_sequence_modify.coord = full_sequence_modify_coord1;
                                index_replace = index_ins+1;
                                ins_count = 0;
                                for i = 1:ins_length
                                    full_sequence_modify.nucleotide(index_replace+ins_count) = split_ins(i,1);
                                    ins_count = ins_count + 1;
                                end
                                
                                row_w = row_w + 1;
                            end
                        end
                        
                        ident_del_sites = strfind(full_sequence_modify.nucleotide,'D');
                        ident_del_sites = find(~cellfun(@isempty,ident_del_sites));
                        full_sequence_modify.nucleotide(ident_del_sites) = [];
                        eval(['haplotype',num2str(n),'.sequence = strjoin(transpose(full_sequence_modify.nucleotide),'''');']);
                    end
                    clearvars full_sequence_modify full_sequence_modify_index1 full_sequence_modify_coord1 full_sequence_modify1
                    if length(reference_seq) > algorithm_type
                        eval(['clearvars haplotype_variants',num2str(n),''])
                    end
                end
            end
            
            counter = 1;
            extra_bases = 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN';
            for h = 1:num_haplotypes
                if isempty(find(num_haplotypes_list==h)) == 0
                    eval(['join_all_together(counter,1) = strcat(cellstr(haplotype',num2str(h),'.sequence),extra_bases);']);
                    counter = counter + 1;
                    if length(reference_seq) > algorithm_type
                        eval(['clearvars haplotype',num2str(h),''])
                    end
                end
            end
            
            full_joined_haplotype_seq = strjoin(join_all_together,'');
            clearvars join_all_together
            index_of_Ns = transpose(strfind(full_joined_haplotype_seq,extra_bases));
            index_of_haplo_start1(1,1) = 1;
            index_of_haplo_start2(1,1) = sgRNA_len;
            index_of_haplo_start1(2:(length(index_of_Ns)+1),1) = index_of_Ns+length(extra_bases);
            index_of_haplo_start2(2:(length(index_of_Ns)+1),1) = index_of_Ns+length(extra_bases)+(sgRNA_len-1);
            
            index_of_haplo_end1 = index_of_haplo_start1-length(extra_bases)-1;
            index_of_haplo_end2 = index_of_haplo_end1-(sgRNA_len-1);
            index_of_haplo_end1(1) = [];
            index_of_haplo_end2(1) = [];
            
            haplo_guide_index1_guides = [];
            haplo_guide_index1_individual = [];
            haplo_guide_index1_guides_only = [];
            haplo_guide_index1_pam_only = [];
            haplo_guide_index1_strand = [];
            haplo_guide_index2_guides = [];
            haplo_guide_index2_individual = [];
            haplo_guide_index2_guides_only = [];
            haplo_guide_index2_pam_only = [];
            haplo_guide_index2_strand = [];
            
            haplo_guide_list_full_temp.guides = [];
            haplo_guide_list_full_temp.bpstart = [];
            haplo_guide_list_full_temp.bpend = [];
            haplo_guide_list_full_temp.strand = [];
            
            disp('Identifying sgRNA for each haplotype...')
            
            if length(PAM_sequence_cell)>1
                for u  = 1:length(PAM_sequence_cell)
                    PAM_sequence_cell_rc{u,1} = seqrcomplement(PAM_sequence_cell{u,1});
                end
            else PAM_sequence_cell_rc = cellstr(seqrcomplement(PAM_sequence_cell{:}));
            end
            
            for m = 1:length(PAM_sequence_cell)
                haplo_guide_index1 = transpose(strfind(full_joined_haplotype_seq,PAM_sequence_cell{m,1}));
            end
            
            if exist('haplo_guide_index1','var') == 1
                if five_prime_PAM == 0
                    dsb_index_remove = haplo_guide_index1-N_seq_in_PAM-4;
                end
                if five_prime_PAM == 1
                    dsb_index_remove = haplo_guide_index1+(N_seq_in_PAM+PAM_len+18-1);
                end
                
                if length(dsb_index_remove)>0
                    counter = 1;
                    for m = 1:length(dsb_index_remove)
                        idx = find((dsb_index_remove(m)-index_of_haplo_end2).*(index_of_haplo_end1-dsb_index_remove(m))>=0, 1);
                        if ~isempty(idx)
                            things_to_remove(m) = 1;
                        else
                            things_to_remove(m) = 0;
                        end
                    end
                    
                    things_to_remove = transpose(things_to_remove);
                    index_to_remove = find(things_to_remove==1);
                    haplo_guide_index1(index_to_remove) = [];
                end
                
                if five_prime_PAM == 0
                    haplo_guide_index1_window_start_temp = haplo_guide_index1-(N_seq_in_PAM+sgRNA_len);
                    haplo_guide_index1_window_end_temp = haplo_guide_index1+(PAM_len-N_seq_in_PAM-1);
                end
                if five_prime_PAM == 1
                    haplo_guide_index1_window_start_temp = haplo_guide_index1;
                    haplo_guide_index1_window_end_temp = haplo_guide_index1+(PAM_len+sgRNA_len-1);
                end
                clearvars haplo_guide_index1
                index_remove_start = find(haplo_guide_index1_window_start_temp <= 0);
                haplo_guide_index1_window_start_temp(index_remove_start) = [];
                haplo_guide_index1_window_end_temp(index_remove_start) = [];
                
                if length(reference_seq) < algorithm_type
                    counter = 0;
                    guide_list = [];
                    for i = 1:(sgRNA_len+PAM_len)
                        eval(['p',num2str(i),' = full_joined_haplotype_seq(haplo_guide_index1_window_start_temp + counter);'])
                        eval(['guide_list = vertcat(guide_list,p',num2str(i),');']);
                        counter = counter + 1;
                    end
                    temp = length(p1);
                    clearvars p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 p16 p17 p18 p19 p20 p21 p22 p23
                elseif length(reference_seq) > algorithm_type
                    counter = 0;
                    guide_list = [];
                    for i = 1:(sgRNA_len+PAM_len)
                        eval(['p',num2str(i),' = full_joined_haplotype_seq(haplo_guide_index1_window_start_temp + counter);'])
                        eval(['guide_list = vertcat(guide_list,p',num2str(i),');']);
                        eval(['temp = length(p',num2str(i),');']);
                        eval(['clearvars p',num2str(i),''])
                        counter = counter + 1;
                    end
                end
                
                for c = 1:temp
                    haplo_guide_list{c,1} = transpose(guide_list(:,c));
                end
            end
            clearvars guide_list haplo_guide_index1_window_start_temp haplo_guide_index1_window_end_temp
            
            for l = 1:length(PAM_sequence_cell)
                haplo_guide_index2 = transpose(strfind(full_joined_haplotype_seq,PAM_sequence_cell_rc{l,1}));
            end
            if exist('haplo_guide_index2','var') == 1
                if five_prime_PAM == 0
                    dsb_index_remove_rc = haplo_guide_index2+(PAM_len+2);
                end
                if five_prime_PAM == 1
                    dsb_index_remove_rc = haplo_guide_index2 + (N_seq_in_PAM+18);
                end
                
                for m = 1:length(dsb_index_remove_rc)
                    idx_rc = find((dsb_index_remove_rc(m)-index_of_haplo_start1).*(index_of_haplo_start2-dsb_index_remove_rc(m))>=0, 1);
                    if ~isempty(idx_rc)
                        things_to_remove_rc(m) = 1;
                    else
                        things_to_remove_rc(m) = 0;
                    end
                end
                
                things_to_remove_rc = transpose(things_to_remove_rc);
                index_to_remove_rc = find(things_to_remove_rc==1);
                haplo_guide_index2(index_to_remove_rc) = [];
                
                if five_prime_PAM == 0
                    haplo_guide_index2_window_start_temp = haplo_guide_index2;
                    haplo_guide_index2_window_end_temp = haplo_guide_index2+(sgRNA_len+PAM_len-1);
                end
                if five_prime_PAM == 1
                    haplo_guide_index2_window_start_temp = haplo_guide_index2-(N_seq_in_PAM+sgRNA_len);
                    haplo_guide_index2_window_end_temp = haplo_guide_index2+(PAM_len-N_seq_in_PAM-1);
                end
                
                clearvars haplo_guide_index2
                index_remove_start_rc = find(haplo_guide_index2_window_start_temp <= 0);
                haplo_guide_index2_window_start_temp(index_remove_start_rc) = [];
                haplo_guide_index2_window_end_temp(index_remove_start_rc) = [];
                
                if length(reference_seq) < algorithm_type
                    counter = 0;
                    guide_list = [];
                    for q = 1:(sgRNA_len+PAM_len)
                        eval(['p',num2str(q),' = full_joined_haplotype_seq(haplo_guide_index2_window_start_temp + counter);'])
                        eval(['guide_list = vertcat(guide_list,p',num2str(q),');']);
                        counter = counter + 1;
                    end
                    temp = length(p1);
                    clearvars p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 p16 p17 p18 p19 p20 p21 p22 p23
                elseif length(reference_seq) > algorithm_type
                    counter = 0;
                    guide_list = [];
                    for q = 1:(sgRNA_len+PAM_len)
                        
                        eval(['p',num2str(q),' = full_joined_haplotype_seq(haplo_guide_index2_window_start_temp + counter);'])
                        eval(['guide_list = vertcat(guide_list,p',num2str(q),');']);
                        eval(['temp = length(p',num2str(q),');']);
                        eval(['clearvars p',num2str(q),''])
                        counter = counter + 1;
                    end
                end
                
                if length(reference_seq) > algorithm_type
                    clearvars full_joined_haplotype_seq
                end
                
                for r = 1:temp
                    haplo_guide_list_rc{r,1} = seqrcomplement(transpose(guide_list(:,r)));
                end
            end
            clearvars guide_list full_joined_haplotype_seq haplo_guide_index2_window_start_temp haplo_guide_index2_window_end_temp
            
            if exist('haplo_guide_list','var') == 1
                haplo_guide_list_strand(1:length(haplo_guide_list),1) = '+';
            end
            if exist('haplo_guide_list_rc','var') == 1
                haplo_guide_list_strand_rc(1:length(haplo_guide_list_rc),1) = '-';
            end
            if exist('haplo_guide_list','var') == 1 & exist('haplo_guide_list_rc','var') == 1
                haplo_guide_list_full_temp.guides = [haplo_guide_list;haplo_guide_list_rc];
                haplo_guide_list_full_temp.strand = [haplo_guide_list_strand;haplo_guide_list_strand_rc];
            elseif exist('haplo_guide_list','var') == 1 & exist('haplo_guide_list_rc','var') == 0
                haplo_guide_list_full_temp.guides = haplo_guide_list;
                haplo_guide_list_full_temp.strand = haplo_guide_list_strand;
            elseif exist('haplo_guide_list','var') == 0 & exist('haplo_guide_list_rc','var') == 1
                haplo_guide_list_full_temp.guides = haplo_guide_list_rc;
                haplo_guide_list_full_temp.strand = haplo_guide_list_strand_rc;
            end
            clearvars haplo_guide_list haplo_guide_list_rc haplo_guide_list_strand haplo_guide_list_strand_rc
            
            disp('Identifying unique sgRNA and counts...')
            
            [haplo_guide_list_full.guides unique_index_haplo] = unique(haplo_guide_list_full_temp.guides);
            haplo_guide_list_full.strand = haplo_guide_list_full_temp.strand(unique_index_haplo);
            index_with_n1 = strfind(haplo_guide_list_full.guides,'N');
            index_with_n = find(~cellfun(@isempty,index_with_n1));
            haplo_guide_list_full.guides(index_with_n) = [];
            haplo_guide_list_full.strand(index_with_n) = [];
            if exist('haplo_guide_list_full')==1 && exist('original_guides_full')==1
                [unique_haplotype_guides_full.guides unique_haplo_index] = setdiff(haplo_guide_list_full.guides,original_guides_full.guides);
                unique_haplotype_guides_full.strand = haplo_guide_list_full.strand(unique_haplo_index);
            elseif exist('haplo_guide_list_full')==1 && exist('original_guides_full')==0
                unique_haplotype_guides_full.guides = haplo_guide_list_full.guides;
                unique_haplotype_guides_full.strand = haplo_guide_list_full.strand;
            end
            if length(unique_haplotype_guides_full.guides)>0
                for k = 1:length(unique_haplotype_guides_full.guides)
                    unique_haplotype_guides_full.count(k,1) = length(find(ismember(haplo_guide_list_full_temp.guides,char(unique_haplotype_guides_full.guides(k,1)))==1));
                end
                unique_haplotype_guides_full.frequency = 100*(unique_haplotype_guides_full.count/num_haplotypes);
                clearvars haplo_guide_list_full_temp
            end
            
            if five_prime_PAM == 0
                for j = 1:length(unique_haplotype_guides_full.guides)
                    unique_haplotype_guides_full.guides_only{j,1} = unique_haplotype_guides_full.guides{j,1}(1:sgRNA_len);
                    unique_haplotype_guides_full.pam_only{j,1} = unique_haplotype_guides_full.guides{j,1}(sgRNA_len+1:length(unique_haplotype_guides_full.guides{j,1}));
                end
            end
            
            if five_prime_PAM == 1
                for j = 1:length(unique_haplotype_guides_full.guides)
                    unique_haplotype_guides_full.guides_only{j,1} = unique_haplotype_guides_full.guides{j,1}((PAM_len+1):(PAM_len+sgRNA_len));
                    unique_haplotype_guides_full.pam_only{j,1} = unique_haplotype_guides_full.guides{j,1}(1:PAM_len);
                end
            end
        end
        
        %% Align Haplotype-Derived Guides
        
        if perform_haplotype_analysis == 1;
            if exist('original_guides_full','var') == 1
                for e = 1:length(unique_haplotype_guides_full.guides)
                    test_sequence = char(unique_haplotype_guides_full.guides(e,1));
                    for g = 1:length(original_guides_full.guides)
                        index_mm(g,1) = nwalign(test_sequence,original_guides_full.guides{g,1});
                    end
                    [min_value(e,1) min_index(e,1)] = max(index_mm);
                end
            end
            
            if length(unique_haplotype_guides_full.strand)>0
                for e = 1:length(unique_haplotype_guides_full.guides_only);
                    if length(strfind(reference_seq,char(unique_haplotype_guides_full.guides_only(e,1))))>0
                        haplotype_exact_match.index(e,1) = min(strfind(reference_seq,char(unique_haplotype_guides_full.guides_only(e,1))));
                        haplotype_exact_match.strand(e,1) = 1;
                    elseif length(strfind(reference_seq,char(seqrcomplement(char(unique_haplotype_guides_full.guides_only(e,1))))))>0
                        haplotype_exact_match.index(e,1) = min(strfind(reference_seq,seqrcomplement(char(unique_haplotype_guides_full.guides_only(e,1)))));
                        haplotype_exact_match.strand(e,1) = 2;
                    elseif length(strfind(reference_seq,char(unique_haplotype_guides_full.guides_only(e,1))))==0
                        haplotype_exact_match.index(e,1) = 0;
                        haplotype_exact_match.strand(e,1) = 0;
                    end
                end
            end
            if exist('haplotype_exact_match','var') == 1
                if length(find(haplotype_exact_match.index>0))>0
                    haplotype_exact_match_list.index = haplotype_exact_match.index(find(haplotype_exact_match.index>0));
                    haplotype_exact_match_list.strand = haplotype_exact_match.strand(find(haplotype_exact_match.index>0));
                end
                
                unique_haplotype_guides_full.pam_create(1:length(haplotype_exact_match.index),1) = {'Variant sgRNA'};
                if length(find(haplotype_exact_match.index>0))>0
                    unique_haplotype_guides_full.pam_create(find(haplotype_exact_match.index>0)) = {'PAM Creation/Alteration'};
                end
                
                if length(find(haplotype_exact_match.index>0))>0
                    row = 1;
                    for n = 1:length(haplotype_exact_match_list.index)
                        if five_prime_PAM == 0 && haplotype_exact_match_list.strand(n,1) == 1 && haplotype_exact_match_list.index(n,1) > 0 && (haplotype_exact_match_list.index(n,1)+sgRNA_len+PAM_len-1)<=length(reference_seq)
                            haplotype_guide_exact_coord_start(row,1) = reference_sequence.coord(haplotype_exact_match_list.index(n,1));
                            haplotype_guide_exact_coord_end(row,1) = haplotype_guide_exact_coord_start(row,1)+sgRNA_len+PAM_len-1;
                            haplotype_guide_exact_index_start(row,1) = haplotype_exact_match_list.index(n,1);
                            haplotype_guide_exact_index_end(row,1) = haplotype_exact_match_list.index(n,1)+sgRNA_len+PAM_len-1;
                            haplotype_guide_exact_dsb_index(row,1) = haplotype_guide_exact_index_start(row,1)+(sgRNA_len-4);
                            haplotype_guide_exact_dsb_coord(row,1) = reference_sequence.coord(haplotype_guide_exact_dsb_index(row,1));
                            haplotype_guide_exact_strand(row,1) = '+';
                        elseif five_prime_PAM == 1 && haplotype_exact_match_list.strand(n,1) == 1 && (haplotype_exact_match_list.index(n,1)-PAM_len)>0 && (haplotype_exact_match_list.index(n,1)+sgRNA_len-1)<=length(reference_seq)
                            haplotype_guide_exact_coord_start(row,1) = reference_sequence.coord(haplotype_exact_match_list.index(n,1))-PAM_len;
                            haplotype_guide_exact_coord_end(row,1) = reference_sequence.coord(haplotype_exact_match_list.index(n,1))+sgRNA_len-1;
                            haplotype_guide_exact_index_start(row,1) = haplotype_exact_match_list.index(n,1)-PAM_len;
                            haplotype_guide_exact_index_end(row,1) = haplotype_exact_match_list.index(n,1)+sgRNA_len-1;
                            haplotype_guide_exact_dsb_index(row,1) = haplotype_guide_exact_index_start(row,1)+(N_seq_in_PAM+PAM_len+18-1);
                            haplotype_guide_exact_dsb_coord(row,1) = reference_sequence.coord(haplotype_guide_exact_dsb_index(row,1));
                            haplotype_guide_exact_strand(row,1) = '+';
                        elseif five_prime_PAM == 0 && haplotype_exact_match_list.strand(n,1) == 2 && (haplotype_exact_match_list.index(n,1)-PAM_len)>0 && (haplotype_exact_match_list.index(n,1)+sgRNA_len-1)<=length(reference_seq)
                            haplotype_guide_exact_coord_start(row,1) = reference_sequence.coord(haplotype_exact_match_list.index(n,1))-PAM_len;
                            haplotype_guide_exact_coord_end(row,1) = haplotype_guide_exact_coord_start(row,1)+sgRNA_len-1;
                            haplotype_guide_exact_index_start(row,1) = haplotype_exact_match_list.index(n,1)-PAM_len;
                            haplotype_guide_exact_index_end(row,1) = haplotype_exact_match_list.index(n,1)+sgRNA_len-1;
                            haplotype_guide_exact_dsb_index(row,1) = haplotype_guide_exact_index_start(row,1)+(PAM_len+2);
                            haplotype_guide_exact_dsb_coord(row,1) = reference_sequence.coord(haplotype_guide_exact_dsb_index(row,1));
                            haplotype_guide_exact_strand(row,1) = '-';
                        elseif five_prime_PAM == 1 && haplotype_exact_match_list.strand(n,1) == 2 && haplotype_exact_match_list.index(n,1) > 0 && (haplotype_exact_match_list.index(n,1)+sgRNA_len+PAM_len-1)<=length(reference_seq)
                            haplotype_guide_exact_coord_start(row,1) = reference_sequence.coord(haplotype_exact_match_list.index(n,1));
                            haplotype_guide_exact_coord_end(row,1) = reference_sequence.coord(haplotype_exact_match_list.index(n,1))+sgRNA_len+PAM_len-1;
                            haplotype_guide_exact_index_start(row,1) = haplotype_exact_match_list.index(n,1);
                            haplotype_guide_exact_index_end(row,1) = haplotype_exact_match_list.index(n,1)+PAM_len+sgRNA_len-1;
                            haplotype_guide_exact_dsb_index(row,1) = haplotype_guide_exact_index_start(row,1)+(N_seq_in_PAM+18);
                            haplotype_guide_exact_dsb_coord(row,1) = reference_sequence.coord(haplotype_guide_exact_dsb_index(row,1));
                            haplotype_guide_exact_strand(row,1) = '-';
                        end
                        row = row + 1;
                    end
                end
                
                if length(find(haplotype_exact_match.index>0))>0;
                    for n = 1:length(haplotype_guide_exact_coord_start)
                        if strmatch(haplotype_guide_exact_strand(n,1),'+') == 1
                            haplotype_exact.guides{n,1} = strjoin(transpose(reference_sequence.nucleotide(haplotype_guide_exact_index_start(n,1):haplotype_guide_exact_index_end(n,1))),'');
                            haplotype_exact.strand{n,1} = '+';
                        elseif strmatch(haplotype_guide_exact_strand(n,1),'-') == 1
                            haplotype_exact.guides{n,1} = seqrcomplement(strjoin(transpose(reference_sequence.nucleotide(haplotype_guide_exact_index_start(n,1):haplotype_guide_exact_index_end(n,1))),''));
                            haplotype_exact.strand{n,1} = '-';
                        end
                        if five_prime_PAM == 0;
                            haplotype_exact.guides_only{n,1} = haplotype_exact.guides{n,1}(1:sgRNA_len);
                            haplotype_exact.pam_only{n,1} = haplotype_exact.guides{n,1}(sgRNA_len+1:length(haplotype_exact.guides{n,1}));
                        elseif five_prime_PAM == 1
                            haplotype_exact.guides_only{n,1} = haplotype_exact.guides{n,1}((PAM_len+1):(PAM_len+sgRNA_len));
                            haplotype_exact.pam_only{n,1} = haplotype_exact.guides{n,1}(1:PAM_len);
                        end
                        haplotype_exact.dsb_index(n,1) = haplotype_guide_exact_dsb_index(n,1);
                        haplotype_exact.dsb_coord(n,1) = haplotype_guide_exact_dsb_coord(n,1);
                        haplotype_exact.guide_index_start(n,1) = haplotype_guide_exact_index_start(n,1);
                        haplotype_exact.guide_index_end(n,1) = haplotype_guide_exact_index_end(n,1);
                        haplotype_exact.guide_coord_start(n,1) = haplotype_guide_exact_coord_start(n,1);
                        haplotype_exact.guide_coord_end(n,1) = haplotype_guide_exact_coord_end(n,1);
                    end
                end
            end
            
            if exist('original_guides_full','var') == 1 && length(unique_haplotype_guides_full.guides)>0
                unique_haplotype_guides_full.n_mm_wt = min_value;
                unique_haplotype_guides_full.guides_wt = original_guides_full.guides(min_index);
                unique_haplotype_guides_full.dsb_coord_wt = original_guides_full.dsb_coord(min_index);
                unique_haplotype_guides_full.guides_only_wt = original_guides_full.guides_only(min_index);
                unique_haplotype_guides_full.pam_only_wt = original_guides_full.pam_only(min_index);
                unique_haplotype_guides_full.strand_wt = original_guides_full.strand(min_index);
                unique_haplotype_guides_full.dsb_index_wt = original_guides_full.dsb_index(min_index);
                unique_haplotype_guides_full.dsb_coord_wt = original_guides_full.dsb_coord(min_index);
                unique_haplotype_guides_full.guide_index_start_wt = original_guides_full.guide_index_start(min_index);
                unique_haplotype_guides_full.guide_index_end_wt = original_guides_full.guide_index_end(min_index);
                unique_haplotype_guides_full.guide_coord_start_wt = original_guides_full.guide_coord_start(min_index);
                unique_haplotype_guides_full.guide_coord_end_wt = original_guides_full.guide_coord_end(min_index);
                if iscell(unique_haplotype_guides_full.n_mm_wt) == 0
                    unique_haplotype_guides_full.n_mm_wt = cellstr(num2str(unique_haplotype_guides_full.n_mm_wt));
                end
            end
            if exist('haplotype_exact_match','var')==1
                unique_haplotype_guides_full.n_mm_wt(find(haplotype_exact_match.index>0)) = {'N/A'};
                if exist('haplotype_exact','var') == 1
                    unique_haplotype_guides_full.guides_wt(find(haplotype_exact_match.index>0)) = haplotype_exact.guides;
                    unique_haplotype_guides_full.dsb_coord_wt(find(haplotype_exact_match.index>0)) = haplotype_exact.dsb_coord;
                    unique_haplotype_guides_full.guides_only_wt(find(haplotype_exact_match.index>0)) = haplotype_exact.guides_only;
                    unique_haplotype_guides_full.pam_only_wt(find(haplotype_exact_match.index>0)) = haplotype_exact.pam_only;
                    unique_haplotype_guides_full.strand_wt(find(haplotype_exact_match.index>0)) = haplotype_exact.strand;
                    unique_haplotype_guides_full.dsb_index_wt(find(haplotype_exact_match.index>0)) = haplotype_exact.dsb_index;
                    unique_haplotype_guides_full.dsb_coord_wt(find(haplotype_exact_match.index>0)) = haplotype_exact.dsb_coord;
                    unique_haplotype_guides_full.guide_index_start_wt(find(haplotype_exact_match.index>0)) = haplotype_exact.guide_index_start;
                    unique_haplotype_guides_full.guide_index_end_wt(find(haplotype_exact_match.index>0)) = haplotype_exact.guide_index_end;
                    unique_haplotype_guides_full.guide_coord_start_wt(find(haplotype_exact_match.index>0)) = haplotype_exact.guide_coord_start;
                    unique_haplotype_guides_full.guide_coord_end_wt(find(haplotype_exact_match.index>0)) = haplotype_exact.guide_coord_end;
                end
            end
            
            if length(unique_haplotype_guides_full.guides)>0
                chr_string = 'chr';
                for m = 1:length(unique_haplotype_guides_full.guides)
                    if strmatch(unique_haplotype_guides_full.strand(m,1),'+') == 1
                        idx = find((unique_haplotype_guides_full.guide_coord_end_wt(m)-full_list_haplo_variants.bpstart).*(full_list_haplo_variants.bpend-unique_haplotype_guides_full.guide_coord_start_wt(m))>=0);
                        num_mult_var_haplo(m,1) = length(idx);
                        num_associated_variants(b,1) = max(num_mult_var_haplo);
                        num_associated_variants_temp = num_associated_variants;
                        num_associated_variants_temp(b,:) = [];
                        num_spots(b,1) = length(num_mult_var_haplo);
                        num_spots_sum(b,1) = sum(num_spots);
                        for e = 1:length(idx);
                            eval(['haplotype_associated_variants',num2str(e),'.bpstart(m,1) = cellstr(num2str(full_list_haplo_variants.bpstart(idx(e))));']);
                            eval(['haplotype_associated_variants',num2str(e),'.bpend(m,1) = cellstr(num2str(full_list_haplo_variants.bpend(idx(e))));']);
                            temp_chr_var = cellstr(strjoin(horzcat(chr_string,cellstr(num2str(full_list_haplo_variants.chr_id(idx(e))))),''));
                            eval(['haplotype_associated_variants',num2str(e),'.chr_id(m,1) = temp_chr_var;']);
                            eval(['haplotype_associated_variants',num2str(e),'.rs_id(m,1) = full_list_haplo_variants.rs_id(idx(e));']);
                            eval(['haplotype_associated_variants',num2str(e),'.ref_snp(m,1) = full_list_haplo_variants.ref_snp(idx(e));']);
                            eval(['haplotype_associated_variants',num2str(e),'.alt_snp(m,1) = full_list_haplo_variants.alt_snp(idx(e));']);
                            eval(['haplotype_associated_variants',num2str(e),'.identifier(m,1) = full_list_haplo_variants.identifier(idx(e));']);
                            eval(['haplotype_associated_variants',num2str(e),'.variant_index_start(m,1) = cellstr(num2str(full_list_haplo_variants.variant_index_start(idx(e))));']);
                            eval(['haplotype_associated_variants',num2str(e),'.variant_index_end(m,1) = cellstr(num2str(full_list_haplo_variants.variant_index_end(idx(e))));']);
                        end
                    elseif strmatch(unique_haplotype_guides_full.strand(m,1),'-') == 1
                        idx = find((unique_haplotype_guides_full.guide_coord_start_wt(m)-full_list_haplo_variants.bpstart).*(full_list_haplo_variants.bpend-unique_haplotype_guides_full.guide_coord_end_wt(m))>=0);
                        num_mult_var_haplo(m,1) = length(idx);
                        num_associated_variants(b,1) = max(num_mult_var_haplo);
                        num_associated_variants_temp = num_associated_variants;
                        num_associated_variants_temp(b,:) = [];
                        
                        num_spots(b,1) = length(num_mult_var_haplo);
                        num_spots_sum(b,1) = sum(num_spots);
                        for e = 1:length(idx);
                            eval(['haplotype_associated_variants',num2str(e),'.bpstart(m,1) = cellstr(num2str(full_list_haplo_variants.bpstart(idx(e))));']);
                            eval(['haplotype_associated_variants',num2str(e),'.bpend(m,1) = cellstr(num2str(full_list_haplo_variants.bpend(idx(e))));']);
                            temp_chr_var = cellstr(strjoin(horzcat(chr_string,cellstr(num2str(full_list_haplo_variants.chr_id(idx(e))))),''));
                            eval(['haplotype_associated_variants',num2str(e),'.chr_id(m,1) = temp_chr_var;']);
                            eval(['haplotype_associated_variants',num2str(e),'.rs_id(m,1) = full_list_haplo_variants.rs_id(idx(e));']);
                            eval(['haplotype_associated_variants',num2str(e),'.ref_snp(m,1) = full_list_haplo_variants.ref_snp(idx(e));']);
                            eval(['haplotype_associated_variants',num2str(e),'.alt_snp(m,1) = full_list_haplo_variants.alt_snp(idx(e));']);
                            eval(['haplotype_associated_variants',num2str(e),'.identifier(m,1) = full_list_haplo_variants.identifier(idx(e));']);
                            eval(['haplotype_associated_variants',num2str(e),'.variant_index_start(m,1) = cellstr(num2str(full_list_haplo_variants.variant_index_start(idx(e))));']);
                            eval(['haplotype_associated_variants',num2str(e),'.variant_index_end(m,1) = cellstr(num2str(full_list_haplo_variants.variant_index_end(idx(e))));']);
                        end
                    end
                end
                
                for e = 1:max(num_mult_var_haplo);
                    if eval(['length(haplotype_associated_variants',num2str(e),'.bpstart) < num_spots(b,1);'])
                        eval(['len_addition = length(haplotype_associated_variants',num2str(e),'.bpstart);']);
                        eval(['haplotype_associated_variants',num2str(e),'.bpstart((len_addition+1):num_spots(b,1),1) = {''''};'])
                        eval(['haplotype_associated_variants',num2str(e),'.bpend((len_addition+1):num_spots(b,1),1) = {''''};'])
                        eval(['haplotype_associated_variants',num2str(e),'.chr_id((len_addition+1):num_spots(b,1),1) = {''''};'])
                        eval(['haplotype_associated_variants',num2str(e),'.rs_id((len_addition+1):num_spots(b,1),1) = {''''};'])
                        eval(['haplotype_associated_variants',num2str(e),'.ref_snp((len_addition+1):num_spots(b,1),1) = {''''};'])
                        eval(['haplotype_associated_variants',num2str(e),'.alt_snp((len_addition+1):num_spots(b,1),1) = {''''};'])
                        eval(['haplotype_associated_variants',num2str(e),'.identifier((len_addition+1):num_spots(b,1),1) = {''''};'])
                        eval(['haplotype_associated_variants',num2str(e),'.variant_index_start((len_addition+1):num_spots(b,1),1) = {''''};'])
                        eval(['haplotype_associated_variants',num2str(e),'.variant_index_end((len_addition+1):num_spots(b,1),1) = {''''};'])
                        
                    end
                end
            end
            
            if twenty_extra_bases == 1 && perform_haplotype_analysis == 1 && exist('original_guides_full','var')==1
                edge_index_left = find(original_guides_full.dsb_coord<(seq_bpstart+sgRNA_len));
                original_guides_full.guides(edge_index_left) = [];
                original_guides_full.guides_only(edge_index_left) = [];
                original_guides_full.pam_only(edge_index_left) = [];
                original_guides_full.strand(edge_index_left) = [];
                original_guides_full.dsb_index(edge_index_left) = [];
                original_guides_full.dsb_coord(edge_index_left) = [];
                original_guides_full.guide_index_start(edge_index_left) = [];
                original_guides_full.guide_index_end(edge_index_left) = [];
                original_guides_full.guide_coord_start(edge_index_left) = [];
                original_guides_full.guide_coord_end(edge_index_left) = [];
                edge_index_right = find(original_guides_full.dsb_coord>(seq_bpend-sgRNA_len));
                original_guides_full.guides(edge_index_right) = [];
                original_guides_full.guides_only(edge_index_right) = [];
                original_guides_full.pam_only(edge_index_right) = [];
                original_guides_full.strand(edge_index_right) = [];
                original_guides_full.dsb_index(edge_index_right) = [];
                original_guides_full.dsb_coord(edge_index_right) = [];
                original_guides_full.guide_index_start(edge_index_right) = [];
                original_guides_full.guide_index_end(edge_index_right) = [];
                original_guides_full.guide_coord_start(edge_index_right) = [];
                original_guides_full.guide_coord_end(edge_index_right) = [];
            end
            
            if perform_haplotype_analysis == 1 && exist('original_guides_full','var') == 1;
                if a == starting_number;
                    for n = 1:length(original_guides_full.guides)
                        zeros_to_add = [];
                        if exist('guide_ID_num_temp','var') == 0
                            guide_ID_num_temp = n;
                        else guide_ID_num_temp = guide_ID_num_temp + 1;
                        end
                        num_zeros_add = length(num2str(guide_ID_num_temp));
                        if num_zeros_add < max_number_guide_ID_digits
                            for j = 1:(max_number_guide_ID_digits-num_zeros_add)
                                zeros_to_add = strcat(zeros_to_add,'0');
                            end
                            original_guides_full.guide_ID_number{n,1} = strcat('sgRNA_',zeros_to_add,num2str(guide_ID_num_temp));
                        elseif num_zeros_add == max_number_guide_ID_digits
                            original_guides_full.guide_ID_number{n,1} = ['sgRNA_',num2str(guide_ID_num_temp),''];
                        end
                    end
                elseif a > starting_number
                    for n = 1:length(original_guides_full.guides)
                        zeros_to_add = [];
                        guide_ID_num_temp = guide_ID_num_temp+1;
                        num_zeros_add = length(num2str(guide_ID_num_temp));
                        if num_zeros_add < max_number_guide_ID_digits
                            for j = 1:(max_number_guide_ID_digits-num_zeros_add)
                                zeros_to_add = strcat(zeros_to_add,'0');
                            end
                            original_guides_full.guide_ID_number{n,1} = strcat('sgRNA_',zeros_to_add,num2str(guide_ID_num_temp));
                        elseif num_zeros_add == max_number_guide_ID_digits
                            original_guides_full.guide_ID_number{n,1} = ['sgRNA_',num2str(guide_ID_num_temp),''];
                        end
                    end
                end
            end
            
            if twenty_extra_bases == 1 && perform_haplotype_analysis == 1 && length(unique_haplotype_guides_full.guides)>0;
                edge_index_left = find(unique_haplotype_guides_full.dsb_coord_wt<(seq_bpstart+sgRNA_len));
                unique_haplotype_guides_full.guides(edge_index_left) = [];
                unique_haplotype_guides_full.guides_only(edge_index_left) = [];
                unique_haplotype_guides_full.pam_only(edge_index_left) = [];
                unique_haplotype_guides_full.strand(edge_index_left) = [];
                unique_haplotype_guides_full.pam_create(edge_index_left) = [];
                unique_haplotype_guides_full.count(edge_index_left) = [];
                unique_haplotype_guides_full.frequency(edge_index_left) = [];
                
                unique_haplotype_guides_full.n_mm_wt(edge_index_left) = [];
                unique_haplotype_guides_full.guides_wt(edge_index_left) = [];
                unique_haplotype_guides_full.guides_only_wt(edge_index_left) = [];
                unique_haplotype_guides_full.pam_only_wt(edge_index_left) = [];
                unique_haplotype_guides_full.strand_wt(edge_index_left) = [];
                unique_haplotype_guides_full.dsb_index_wt(edge_index_left) = [];
                unique_haplotype_guides_full.dsb_coord_wt(edge_index_left) = [];
                unique_haplotype_guides_full.guide_index_start_wt(edge_index_left) = [];
                unique_haplotype_guides_full.guide_index_end_wt(edge_index_left) = [];
                unique_haplotype_guides_full.guide_coord_start_wt(edge_index_left) = [];
                unique_haplotype_guides_full.guide_coord_end_wt(edge_index_left) = [];
                
                for e = 1:num_associated_variants(b,1)
                    eval(['haplotype_associated_variants',num2str(e),'.bpstart(edge_index_left) = [];']);
                    eval(['haplotype_associated_variants',num2str(e),'.bpend(edge_index_left) = [];']);
                    eval(['haplotype_associated_variants',num2str(e),'.rs_id(edge_index_left) = [];']);
                    eval(['haplotype_associated_variants',num2str(e),'.ref_snp(edge_index_left) = [];']);
                    eval(['haplotype_associated_variants',num2str(e),'.alt_snp(edge_index_left) = [];']);
                    eval(['haplotype_associated_variants',num2str(e),'.identifier(edge_index_left) = [];']);
                    eval(['haplotype_associated_variants',num2str(e),'.variant_index_start(edge_index_left) = [];']);
                    eval(['haplotype_associated_variants',num2str(e),'.variant_index_end(edge_index_left) = [];']);
                    eval(['haplotype_associated_variants',num2str(e),'.chr_id(edge_index_left) = [];']);
                end
                
                edge_index_right = find(unique_haplotype_guides_full.dsb_coord_wt>(seq_bpend-sgRNA_len));
                unique_haplotype_guides_full.guides(edge_index_right) = [];
                unique_haplotype_guides_full.guides_only(edge_index_right) = [];
                unique_haplotype_guides_full.pam_only(edge_index_right) = [];
                unique_haplotype_guides_full.strand(edge_index_right) = [];
                unique_haplotype_guides_full.pam_create(edge_index_right) = [];
                unique_haplotype_guides_full.count(edge_index_right) = [];
                unique_haplotype_guides_full.frequency(edge_index_right) = [];
                
                unique_haplotype_guides_full.n_mm_wt(edge_index_right) = [];
                unique_haplotype_guides_full.guides_wt(edge_index_right) = [];
                unique_haplotype_guides_full.guides_only_wt(edge_index_right) = [];
                unique_haplotype_guides_full.pam_only_wt(edge_index_right) = [];
                unique_haplotype_guides_full.strand_wt(edge_index_right) = [];
                unique_haplotype_guides_full.dsb_index_wt(edge_index_right) = [];
                unique_haplotype_guides_full.dsb_coord_wt(edge_index_right) = [];
                unique_haplotype_guides_full.guide_index_start_wt(edge_index_right) = [];
                unique_haplotype_guides_full.guide_index_end_wt(edge_index_right) = [];
                unique_haplotype_guides_full.guide_coord_start_wt(edge_index_right) = [];
                unique_haplotype_guides_full.guide_coord_end_wt(edge_index_right) = [];
                
                num_spots(b,1) = num_spots(b,1)-length(edge_index_right)-length(edge_index_left);
                num_spots_sum(b,1) = num_spots_sum(b,1)-length(edge_index_right)-length(edge_index_left);
                num_associated_variants_batch(a,1) = max(num_associated_variants);
                
                for e = 1:num_associated_variants(b,1)
                    eval(['haplotype_associated_variants',num2str(e),'.bpstart(edge_index_right) = [];']);
                    eval(['haplotype_associated_variants',num2str(e),'.bpend(edge_index_right) = [];']);
                    eval(['haplotype_associated_variants',num2str(e),'.rs_id(edge_index_right) = [];']);
                    eval(['haplotype_associated_variants',num2str(e),'.ref_snp(edge_index_right) = [];']);
                    eval(['haplotype_associated_variants',num2str(e),'.alt_snp(edge_index_right) = [];']);
                    eval(['haplotype_associated_variants',num2str(e),'.identifier(edge_index_right) = [];']);
                    eval(['haplotype_associated_variants',num2str(e),'.variant_index_start(edge_index_right) = [];']);
                    eval(['haplotype_associated_variants',num2str(e),'.variant_index_end(edge_index_right) = [];']);
                    eval(['haplotype_associated_variants',num2str(e),'.chr_id(edge_index_right) = [];']);
                end
            end
            
            if twenty_extra_bases == 0 && perform_haplotype_analysis == 1 && length(unique_haplotype_guides_full.guides)>0;
                num_associated_variants_batch(a,1) = max(num_associated_variants);
            end
            
            if length(unique_haplotype_guides_full.guides)>0
                for e = 1:length(unique_haplotype_guides_full.guides_only)
                    if length(strfind(reference_seq,char(unique_haplotype_guides_full.guides_only(e,1))))>0
                        unique_haplotype_guides_full.guide_ID_number_wt(e,1) = {'N/A'};
                    elseif length(strfind(reference_seq,seqrcomplement(char(unique_haplotype_guides_full.guides_only(e,1)))))>0
                        unique_haplotype_guides_full.guide_ID_number_wt(e,1) = {'N/A'};
                    elseif length(find(ismember(original_guides_full.guides_only,unique_haplotype_guides_full.guides_only_wt(e,1))>0))>0
                        unique_haplotype_guides_full.guide_ID_number_wt(e,1) = original_guides_full.guide_ID_number(min(find(ismember(original_guides_full.guides_only,unique_haplotype_guides_full.guides_only_wt(e,1))>0)));
                    elseif length(find(ismember(original_guides_full.guides_only,unique_haplotype_guides_full.guides_only_wt(e,1))>0))>0
                        unique_haplotype_guides_full.guide_ID_number_wt(e,1) = original_guides_full.guide_ID_number(min(find(ismember(original_guides_full.guides_only,seqrcomplement(char(unique_haplotype_guides_full.guides_only_wt(e,1))))>0)));
                    end
                end
                
                if exist('guide_ID_num_temp','var') == 0
                    guide_ID_num_temp = 0;
                end
                for n = 1:length(unique_haplotype_guides_full.guides)
                    zeros_to_add = [];
                    guide_ID_num_temp = guide_ID_num_temp+1;
                    num_zeros_add = length(num2str(guide_ID_num_temp));
                    if num_zeros_add < max_number_guide_ID_digits
                        for j = 1:(max_number_guide_ID_digits-num_zeros_add)
                            zeros_to_add = strcat(zeros_to_add,'0');
                        end
                        unique_haplotype_guides_full.guide_ID_number{n,1} = strcat('sgRNA_',zeros_to_add,num2str(guide_ID_num_temp));
                    elseif num_zeros_add == max_number_guide_ID_digits
                        unique_haplotype_guides_full.guide_ID_number{n,1} = ['sgRNA_',num2str(guide_ID_num_temp),''];
                    end
                end
                
            end
        end
        
        %% Variant Analysis
        
        if perform_variant_analysis == 1 || perform_wgs_analysis == 1
            disp('Variant analysis...')
            
            var_guide_index1 = [];
            var_guide_index2 = [];
            variant_analysis.guides = [];
            variant_analysis.ref_snp = [];
            variant_analysis.alt_snp = [];
            variant_analysis.ref_snp = [];
            variant_analysis.alt_snp = [];
            variant_analysis.from_seq = [];
            variant_analysis.variant_index_start = [];
            variant_analysis.variant_index_end = [];
            variant_analysis.variant_coord_start = [];
            variant_analysis.variant_coord_end = [];
            variant_analysis.dsb_coord = [];
            variant_analysis.dsb_index = [];
            variant_analysis.gmaf = [];
            variant_analysis.rs_id = [];
            variant_analysis.identifier = [];
            variant_analysis.variant_len = [];
            variant_analysis.type = [];
            variant_analysis_rc.guides = [];
            variant_analysis_rc.ref_snp = [];
            variant_analysis_rc.alt_snp = [];
            variant_analysis_rc.ref_snp = [];
            variant_analysis_rc.alt_snp = [];
            variant_analysis_rc.from_seq = [];
            variant_analysis_rc.variant_index_start = [];
            variant_analysis_rc.variant_index_end = [];
            variant_analysis_rc.variant_coord_start = [];
            variant_analysis_rc.variant_coord_end = [];
            variant_analysis_rc.dsb_coord = [];
            variant_analysis_rc.dsb_index = [];
            variant_analysis_rc.gmaf = [];
            variant_analysis_rc.rs_id = [];
            variant_analysis_rc.identifier = [];
            variant_analysis_rc.variant_len = [];
            variant_analysis_rc.type = [];
            
            if perform_variant_analysis == 1 && perform_wgs_analysis == 0
                for n = 1:num_windows
                    eval(['var_iterations = length(window_sequence_var',num2str(n),'.window);']);
                    if var_iterations > 1
                        for i = 2:var_iterations
                            eval(['variant_window_seq = window_sequence_var',num2str(n),'.window(i,1);']);
                            eval(['variant_window_index_start = window_sequence_var',num2str(n),'.window_start;']);
                            eval(['variant_window_index_end = window_sequence_var',num2str(n),'.window_end;']);
                            eval(['variant_window_ref_snp = window_sequence_var',num2str(n),'.ref_snp(i,1);']);
                            eval(['variant_window_alt_snp = window_sequence_var',num2str(n),'.alt_snp(i,1);']);
                            eval(['variant_window_from_seq = window_sequence_var',num2str(n),'.from_seq(i,1);']);
                            eval(['variant_window_variant_index_start = window_sequence_var',num2str(n),'.variant_index_start(i,1);']);
                            eval(['variant_window_variant_index_end = window_sequence_var',num2str(n),'.variant_index_end(i,1);']);
                            eval(['variant_window_variant_coord_start = window_sequence_var',num2str(n),'.variant_coord_start(i,1);']);
                            eval(['variant_window_variant_coord_end = window_sequence_var',num2str(n),'.variant_coord_end(i,1);']);
                            eval(['variant_window_rs_id = window_sequence_var',num2str(n),'.rs_id(i,1);']);
                            eval(['variant_window_identifier = window_sequence_var',num2str(n),'.identifier(i,1);']);
                            eval(['variant_window_variant_len = window_sequence_var',num2str(n),'.variant_len(i,1);']);
                            eval(['variant_window_type = window_sequence_var',num2str(n),'.type(i,1);']);
                            variant_window_seq = char(variant_window_seq);
                            variant_window_seq_index = reference_sequence.number(variant_window_index_start:variant_window_index_end);
                            var_guide_index1 = [];
                            var_guide_index2 = [];
                            for m = 1:length(PAM_sequence_cell)
                                var_guide_index1_temp_index = transpose(strfind(variant_window_seq,PAM_sequence_cell{m,1}));
                                var_guide_index1 = [var_guide_index1;var_guide_index1_temp_index];
                            end
                            if five_prime_PAM == 0
                                var_guide_index1_start_temp = var_guide_index1-(N_seq_in_PAM+sgRNA_len);
                                var_guide_index1_end_temp = var_guide_index1+(PAM_len-N_seq_in_PAM-1);
                            elseif five_prime_PAM == 1
                                var_guide_index1_start_temp = var_guide_index1;
                                var_guide_index1_end_temp = var_guide_index1+(PAM_len+sgRNA_len-1);
                            end
                            for k = 1:length(var_guide_index1)
                                if var_guide_index1_start_temp(k,1) > 0 && var_guide_index1_end_temp(k,1)<=length(variant_window_seq)
                                    variant_analysis.guides = [variant_analysis.guides;variant_window_seq(var_guide_index1_start_temp(k,1):var_guide_index1_end_temp(k,1))];
                                    variant_analysis.ref_snp = [variant_analysis.ref_snp;variant_window_ref_snp];
                                    variant_analysis.alt_snp = [variant_analysis.alt_snp;variant_window_alt_snp];
                                    variant_analysis.from_seq = [variant_analysis.from_seq;variant_window_from_seq];
                                    variant_analysis.variant_index_start = [variant_analysis.variant_index_start;variant_window_variant_index_start];
                                    variant_analysis.variant_index_end = [variant_analysis.variant_index_end;variant_window_variant_index_end];
                                    variant_analysis.variant_coord_start = [variant_analysis.variant_coord_start;variant_window_variant_coord_start];
                                    variant_analysis.variant_coord_end = [variant_analysis.variant_coord_end;variant_window_variant_coord_end];
                                    variant_analysis.rs_id = [variant_analysis.rs_id;variant_window_rs_id];
                                    variant_analysis.identifier = [variant_analysis.identifier;variant_window_identifier];
                                    variant_analysis.variant_len = [variant_analysis.variant_len;variant_window_variant_len];
                                    variant_analysis.type = [variant_analysis.type;variant_window_type];
                                end
                            end
                            if exist('variant_guide_list','var') == 1
                                var_guide_index1_guides = [var_guide_index1_guides;variant_guide_list];
                            end
                            for l = 1:length(PAM_sequence_cell)
                                var_guide_index2_temp_index = transpose(strfind(variant_window_seq,PAM_sequence_cell_rc{l,1}));
                                var_guide_index2 = [var_guide_index2;var_guide_index2_temp_index];
                            end
                            if five_prime_PAM == 0
                                var_guide_index2_end_temp = var_guide_index2;
                                var_guide_index2_start_temp = var_guide_index2+(sgRNA_len+PAM_len-1);
                            elseif five_prime_PAM == 1
                                var_guide_index2_end_temp = var_guide_index2-(N_seq_in_PAM+sgRNA_len);
                                var_guide_index2_start_temp = var_guide_index2+PAM_len-N_seq_in_PAM-1;
                            end
                            
                            for w = 1:length(var_guide_index2)
                                if var_guide_index2_end_temp(w,1) > 0 && var_guide_index2_start_temp(w,1)<=length(variant_window_seq)
                                    variant_analysis_rc.guides = [variant_analysis_rc.guides;seqrcomplement(variant_window_seq(var_guide_index2_end_temp(w,1):var_guide_index2_start_temp(w,1)))];
                                    variant_analysis_rc.ref_snp = [variant_analysis_rc.ref_snp;variant_window_ref_snp];
                                    variant_analysis_rc.alt_snp = [variant_analysis_rc.alt_snp;variant_window_alt_snp];
                                    variant_analysis_rc.from_seq = [variant_analysis_rc.from_seq;variant_window_from_seq];
                                    variant_analysis_rc.variant_index_start = [variant_analysis_rc.variant_index_start;variant_window_variant_index_start];
                                    variant_analysis_rc.variant_index_end = [variant_analysis_rc.variant_index_end;variant_window_variant_index_end];
                                    variant_analysis_rc.variant_coord_start = [variant_analysis_rc.variant_coord_start;variant_window_variant_coord_start];
                                    variant_analysis_rc.variant_coord_end = [variant_analysis_rc.variant_coord_end;variant_window_variant_coord_end];
                                    variant_analysis_rc.rs_id = [variant_analysis_rc.rs_id;variant_window_rs_id];
                                    variant_analysis_rc.identifier = [variant_analysis_rc.identifier;variant_window_identifier];
                                    variant_analysis_rc.variant_len = [variant_analysis_rc.variant_len;variant_window_variant_len];
                                    variant_analysis_rc.type = [variant_analysis_rc.type;variant_window_type];
                                end
                            end
                            if exist('variant_guide_list_rc','var') == 1
                                var_guide_index2_guides = [var_guide_index2_guides;variant_guide_list_rc];
                            end
                            clearvars var_guide_index1_temp_index var_guide_index1_temp var_guide_index1 variant_window_seq variant_window_seq_rc
                            clearvars var_guide_index2_temp_index var_guide_index2_temp var_guide_index2 variant_window_seq variant_window_seq_rc
                            clearvars variant_window_ref_snp variant_window_alt_snp variant_window_from_seq variant_window_variant_index_start
                            clearvars variant_window_variant_index_end variant_window_variant_coord_start variant_window_variant_coord_end
                            clearvars variant_window_gmaf variant_window_rs_id variant_window_identifier variant_window_variant_len variant_window_type
                        end
                        clearvars var_iterations variant_window_seq_index variant_guide_list variant_guide_list_rc
                    end
                end
            end
            
            if perform_variant_analysis == 0 && perform_wgs_analysis == 1
                for n = 1:num_windows
                    eval(['var_iterations = length(window_sequence_var',num2str(n),'.window);']);
                    if var_iterations > 1
                        for i = var_iterations:var_iterations
                            eval(['variant_window_seq = window_sequence_var',num2str(n),'.window(i,1);']);
                            eval(['variant_window_index_start = window_sequence_var',num2str(n),'.window_start;']);
                            eval(['variant_window_index_end = window_sequence_var',num2str(n),'.window_end;']);
                            eval(['variant_window_ref_snp = window_sequence_var',num2str(n),'.ref_snp(i,1);']);
                            eval(['variant_window_alt_snp = window_sequence_var',num2str(n),'.alt_snp(i,1);']);
                            eval(['variant_window_from_seq = window_sequence_var',num2str(n),'.from_seq(i,1);']);
                            eval(['variant_window_variant_index_start = window_sequence_var',num2str(n),'.variant_index_start(i,1);']);
                            eval(['variant_window_variant_index_end = window_sequence_var',num2str(n),'.variant_index_end(i,1);']);
                            eval(['variant_window_variant_coord_start = window_sequence_var',num2str(n),'.variant_coord_start(i,1);']);
                            eval(['variant_window_variant_coord_end = window_sequence_var',num2str(n),'.variant_coord_end(i,1);']);
                            eval(['variant_window_rs_id = window_sequence_var',num2str(n),'.rs_id(i,1);']);
                            eval(['variant_window_identifier = window_sequence_var',num2str(n),'.identifier(i,1);']);
                            eval(['variant_window_variant_len = window_sequence_var',num2str(n),'.variant_len(i,1);']);
                            eval(['variant_window_type = window_sequence_var',num2str(n),'.type(i,1);']);
                            variant_window_seq = char(variant_window_seq);
                            variant_window_seq_index = reference_sequence.number(variant_window_index_start:variant_window_index_end);
                            var_guide_index1 = [];
                            var_guide_index2 = [];
                            for m = 1:length(PAM_sequence_cell)
                                var_guide_index1_temp_index = transpose(strfind(variant_window_seq,PAM_sequence_cell{m,1}));
                                var_guide_index1 = [var_guide_index1;var_guide_index1_temp_index];
                            end
                            if five_prime_PAM == 0
                                var_guide_index1_start_temp = var_guide_index1-(N_seq_in_PAM+sgRNA_len);
                                var_guide_index1_end_temp = var_guide_index1+(PAM_len-N_seq_in_PAM-1);
                            elseif five_prime_PAM == 1
                                var_guide_index1_start_temp = var_guide_index1;
                                var_guide_index1_end_temp = var_guide_index1+(PAM_len+sgRNA_len-1);
                            end
                            for k = 1:length(var_guide_index1)
                                if var_guide_index1_start_temp(k,1) > 0 && var_guide_index1_end_temp(k,1)<=length(variant_window_seq)
                                    variant_analysis.guides = [variant_analysis.guides;variant_window_seq(var_guide_index1_start_temp(k,1):var_guide_index1_end_temp(k,1))];
                                    variant_analysis.ref_snp = [variant_analysis.ref_snp;variant_window_ref_snp];
                                    variant_analysis.alt_snp = [variant_analysis.alt_snp;variant_window_alt_snp];
                                    variant_analysis.from_seq = [variant_analysis.from_seq;variant_window_from_seq];
                                    variant_analysis.variant_index_start = [variant_analysis.variant_index_start;variant_window_variant_index_start];
                                    variant_analysis.variant_index_end = [variant_analysis.variant_index_end;variant_window_variant_index_end];
                                    variant_analysis.variant_coord_start = [variant_analysis.variant_coord_start;variant_window_variant_coord_start];
                                    variant_analysis.variant_coord_end = [variant_analysis.variant_coord_end;variant_window_variant_coord_end];
                                    variant_analysis.rs_id = [variant_analysis.rs_id;variant_window_rs_id];
                                    variant_analysis.identifier = [variant_analysis.identifier;variant_window_identifier];
                                    variant_analysis.variant_len = [variant_analysis.variant_len;variant_window_variant_len];
                                    variant_analysis.type = [variant_analysis.type;variant_window_type];
                                end
                            end
                            if exist('variant_guide_list','var') == 1
                                var_guide_index1_guides = [var_guide_index1_guides;variant_guide_list];
                            end
                            for l = 1:length(PAM_sequence_cell)
                                var_guide_index2_temp_index = transpose(strfind(variant_window_seq,PAM_sequence_cell_rc{l,1}));
                                var_guide_index2 = [var_guide_index2;var_guide_index2_temp_index];
                            end
                            if five_prime_PAM == 0
                                var_guide_index2_end_temp = var_guide_index2;
                                var_guide_index2_start_temp = var_guide_index2+(sgRNA_len+PAM_len-1);
                            elseif five_prime_PAM == 1
                                var_guide_index2_end_temp = var_guide_index2-(N_seq_in_PAM+sgRNA_len);
                                var_guide_index2_start_temp = var_guide_index2+PAM_len-N_seq_in_PAM-1;
                            end
                            
                            for w = 1:length(var_guide_index2)
                                if var_guide_index2_end_temp(w,1) > 0 && var_guide_index2_start_temp(w,1)<=length(variant_window_seq)
                                    variant_analysis_rc.guides = [variant_analysis_rc.guides;seqrcomplement(variant_window_seq(var_guide_index2_end_temp(w,1):var_guide_index2_start_temp(w,1)))];
                                    variant_analysis_rc.ref_snp = [variant_analysis_rc.ref_snp;variant_window_ref_snp];
                                    variant_analysis_rc.alt_snp = [variant_analysis_rc.alt_snp;variant_window_alt_snp];
                                    variant_analysis_rc.from_seq = [variant_analysis_rc.from_seq;variant_window_from_seq];
                                    variant_analysis_rc.variant_index_start = [variant_analysis_rc.variant_index_start;variant_window_variant_index_start];
                                    variant_analysis_rc.variant_index_end = [variant_analysis_rc.variant_index_end;variant_window_variant_index_end];
                                    variant_analysis_rc.variant_coord_start = [variant_analysis_rc.variant_coord_start;variant_window_variant_coord_start];
                                    variant_analysis_rc.variant_coord_end = [variant_analysis_rc.variant_coord_end;variant_window_variant_coord_end];
                                    variant_analysis_rc.rs_id = [variant_analysis_rc.rs_id;variant_window_rs_id];
                                    variant_analysis_rc.identifier = [variant_analysis_rc.identifier;variant_window_identifier];
                                    variant_analysis_rc.variant_len = [variant_analysis_rc.variant_len;variant_window_variant_len];
                                    variant_analysis_rc.type = [variant_analysis_rc.type;variant_window_type];
                                end
                            end
                            if exist('variant_guide_list_rc','var') == 1
                                var_guide_index2_guides = [var_guide_index2_guides;variant_guide_list_rc];
                            end
                            clearvars var_guide_index1_temp_index var_guide_index1_temp var_guide_index1 variant_window_seq variant_window_seq_rc
                            clearvars var_guide_index2_temp_index var_guide_index2_temp var_guide_index2 variant_window_seq variant_window_seq_rc
                            clearvars variant_window_ref_snp variant_window_alt_snp variant_window_from_seq variant_window_variant_index_start
                            clearvars variant_window_variant_index_end variant_window_variant_coord_start variant_window_variant_coord_end
                            clearvars variant_window_gmaf variant_window_rs_id variant_window_identifier variant_window_variant_len variant_window_type
                        end
                        clearvars var_iterations variant_window_seq_index variant_guide_list variant_guide_list_rc
                    end
                end
            end
            
            
            if length(variant_analysis.variant_index_start)>0
                if length(variant_analysis.guides) > 0
                    variant_analysis.guides = cellstr(variant_analysis.guides);
                    for n = 1:length(variant_analysis.guides)
                        if five_prime_PAM == 0
                            variant_analysis.guides_only{n,1} = variant_analysis.guides{n,1}(1:sgRNA_len);
                            variant_analysis.pam_only{n,1} = variant_analysis.guides{n,1}(sgRNA_len+1:length(variant_analysis.guides{n,1}));
                        elseif five_prime_PAM == 1
                            variant_analysis.guides_only{n,1} = variant_analysis.guides{n,1}((PAM_len+1):(PAM_len+sgRNA_len));
                            variant_analysis.pam_only{n,1} = variant_analysis.guides{n,1}(1:PAM_len);
                        end
                        variant_analysis.strand{n,1} = '+';
                    end
                elseif length(variant_analysis.guides) == 0
                    variant_analysis.guides_only = [];
                    variant_analysis.pam_only = [];
                    variant_analysis.strand = [];
                end
                if length(variant_analysis_rc.guides) > 0
                    variant_analysis_rc.guides = cellstr(variant_analysis_rc.guides);
                    for n = 1:length(variant_analysis_rc.guides)
                        if five_prime_PAM == 0
                            variant_analysis_rc.guides_only{n,1} = variant_analysis_rc.guides{n,1}(1:sgRNA_len);
                            variant_analysis_rc.pam_only{n,1} = variant_analysis_rc.guides{n,1}(sgRNA_len+1:length(variant_analysis_rc.guides{n,1}));
                        elseif five_prime_PAM == 1
                            variant_analysis_rc.guides_only{n,1} = variant_analysis_rc.guides{n,1}((PAM_len+1):(PAM_len+sgRNA_len));
                            variant_analysis_rc.pam_only{n,1} = variant_analysis_rc.guides{n,1}(1:PAM_len);
                        end
                        variant_analysis_rc.strand{n,1} = '-';
                    end
                elseif length(variant_analysis_rc.guides) == 0
                    variant_analysis_rc.guides_only = [];
                    variant_analysis_rc.pam_only = [];
                    variant_analysis_rc.strand = [];
                end
                variant_guides_full_temp.guides = [variant_analysis.guides;variant_analysis_rc.guides];
                variant_guides_full_temp.guides_only = [variant_analysis.guides_only;variant_analysis_rc.guides_only];
                variant_guides_full_temp.pam_only = [variant_analysis.pam_only;variant_analysis_rc.pam_only];
                variant_guides_full_temp.strand = [variant_analysis.strand;variant_analysis_rc.strand];
                variant_guides_full_temp.variant_index_start = [variant_analysis.variant_index_start;variant_analysis_rc.variant_index_start];
                variant_guides_full_temp.variant_index_end = [variant_analysis.variant_index_end;variant_analysis_rc.variant_index_end];
                variant_guides_full_temp.variant_coord_start = [variant_analysis.variant_coord_start;variant_analysis_rc.variant_coord_start];
                variant_guides_full_temp.variant_coord_end = [variant_analysis.variant_coord_end;variant_analysis_rc.variant_coord_end];
                variant_guides_full_temp.ref_snp = [variant_analysis.ref_snp;variant_analysis_rc.ref_snp];
                variant_guides_full_temp.alt_snp = [variant_analysis.alt_snp;variant_analysis_rc.alt_snp];
                variant_guides_full_temp.from_seq = [variant_analysis.from_seq;variant_analysis_rc.from_seq];
                variant_guides_full_temp.rs_id = [variant_analysis.rs_id;variant_analysis_rc.rs_id];
                variant_guides_full_temp.identifier = [variant_analysis.identifier;variant_analysis_rc.identifier];
                variant_guides_full_temp.type = [variant_analysis.type;variant_analysis_rc.type];
                variant_guides_full_temp.variant_len = [variant_analysis.variant_len;variant_analysis_rc.variant_len];
                
                clearvars unique_index_full
                [variant_guides_full.guides unique_index_full] = unique(variant_guides_full_temp.guides);
                variant_guides_full.guides_only = variant_guides_full_temp.guides_only(unique_index_full);
                variant_guides_full.pam_only = variant_guides_full_temp.pam_only(unique_index_full);
                variant_guides_full.strand = variant_guides_full_temp.strand(unique_index_full);
                variant_guides_full.variant_index_start = variant_guides_full_temp.variant_index_start(unique_index_full);
                variant_guides_full.variant_index_end = variant_guides_full_temp.variant_index_end(unique_index_full);
                variant_guides_full.variant_coord_start = variant_guides_full_temp.variant_coord_start(unique_index_full);
                variant_guides_full.variant_coord_end = variant_guides_full_temp.variant_coord_end(unique_index_full);
                variant_guides_full.ref_snp = variant_guides_full_temp.ref_snp(unique_index_full);
                variant_guides_full.alt_snp = variant_guides_full_temp.alt_snp(unique_index_full);
                variant_guides_full.from_seq = variant_guides_full_temp.from_seq(unique_index_full);
                variant_guides_full.rs_id = variant_guides_full_temp.rs_id(unique_index_full);
                variant_guides_full.identifier = variant_guides_full_temp.identifier(unique_index_full);
                variant_guides_full.type = variant_guides_full_temp.type(unique_index_full);
                variant_guides_full.variant_len = variant_guides_full_temp.variant_len(unique_index_full);
                
                [unique_variant_guides_full.guides unique_var_index] = setdiff(variant_guides_full.guides,original_guides_full.guides);
                unique_variant_guides_full.guides_only = variant_guides_full.guides_only(unique_var_index);
                unique_variant_guides_full.pam_only = variant_guides_full.pam_only(unique_var_index);
                unique_variant_guides_full.strand = variant_guides_full.strand(unique_var_index);
                unique_variant_guides_full.from_seq = variant_guides_full.from_seq(unique_var_index);
                unique_variant_guides_full.identifier = variant_guides_full.identifier(unique_var_index);
                unique_variant_guides_full.type = variant_guides_full.type(unique_var_index);
                unique_variant_guides_full.variant_len = variant_guides_full.variant_len(unique_var_index);
                
                if perform_variant_analysis == 1 && length(unique_variant_guides_full.guides)>0 || perform_wgs_analysis == 1 && length(unique_variant_guides_full.guides)>0;
                    for e = 1:length(unique_variant_guides_full.guides)
                        test_sequence = char(unique_variant_guides_full.guides(e,1));
                        if length(original_guides_full.guides) >0
                        for g = 1:length(original_guides_full.guides)
                            index_mm(g,1) = nwalign(test_sequence,original_guides_full.guides{g,1});
                        end
                        [min_value(e,1) min_index(e,1)] = max(index_mm);
                        elseif length(original_guides_full.guides) == 0
                            min_value(e,1) = 0;
                            min_index(e,1) = 0;
                        end
                    end
                    
                    for e = 1:length(unique_variant_guides_full.guides_only);
                        if length(strfind(reference_seq,char(unique_variant_guides_full.guides_only(e,1))))>0
                            variant_exact_match.index(e,1) = min(strfind(reference_seq,char(unique_variant_guides_full.guides_only(e,1))));
                            variant_exact_match.strand(e,1) = 1;
                        elseif length(strfind(reference_seq,char(seqrcomplement(char(unique_variant_guides_full.guides_only(e,1))))))>0
                            variant_exact_match.index(e,1) = min(strfind(reference_seq,seqrcomplement(char(unique_variant_guides_full.guides_only(e,1)))));
                            variant_exact_match.strand(e,1) = 2;
                        elseif length(strfind(reference_seq,char(unique_variant_guides_full.guides_only(e,1))))==0
                            variant_exact_match.index(e,1) = 0;
                            variant_exact_match.strand(e,1) = 0;
                        end
                    end
                    
                    
                    if length(find(variant_exact_match.index>0))>0
                        variant_exact_match_list.index = variant_exact_match.index(find(variant_exact_match.index>0));
                        variant_exact_match_list.strand = variant_exact_match.strand(find(variant_exact_match.index>0));
                    end
                    
                    unique_variant_guides_full.pam_create(1:length(variant_exact_match.index),1) = {'PAM Creation/Alteration'};
                    if length(find(variant_exact_match.index==0))>0
                        unique_variant_guides_full.pam_create(find(variant_exact_match.index==0)) = {'Variant sgRNA'};
                    end
                    
                    if length(find(variant_exact_match.index>0))>0
                        row = 1;
                        for n = 1:length(variant_exact_match_list.index)
                            if five_prime_PAM == 0 && variant_exact_match_list.strand(n,1) == 1 && variant_exact_match_list.index(n,1) > 0
                                variant_guide_exact_coord_start(row,1) = reference_sequence.coord(variant_exact_match_list.index(n,1));
                                variant_guide_exact_coord_end(row,1) = variant_guide_exact_coord_start(row,1)+sgRNA_len+PAM_len-1;
                                variant_guide_exact_index_start(row,1) = variant_exact_match_list.index(n,1);
                                variant_guide_exact_index_end(row,1) = variant_exact_match_list.index(n,1)+sgRNA_len+PAM_len-1;
                                variant_guide_exact_dsb_index(row,1) = variant_guide_exact_index_start(row,1)+(sgRNA_len-4);
                                variant_guide_exact_dsb_coord(row,1) = reference_sequence.coord(variant_guide_exact_dsb_index(row,1));
                                variant_guide_exact_strand(row,1) = '+';
                            elseif five_prime_PAM == 1 && variant_exact_match_list.strand(n,1) == 1 && (variant_exact_match_list.index(n,1)-PAM_len)>0
                                variant_guide_exact_coord_start(row,1) = reference_sequence.coord(variant_exact_match_list.index(n,1))-PAM_len;
                                variant_guide_exact_coord_end(row,1) = reference_sequence.coord(variant_exact_match_list.index(n,1))+sgRNA_len-1;
                                variant_guide_exact_index_start(row,1) = variant_exact_match_list.index(n,1)-PAM_len;
                                variant_guide_exact_index_end(row,1) = variant_exact_match_list.index(n,1)+sgRNA_len-1;
                                variant_guide_exact_dsb_index(row,1) = variant_guide_exact_index_start(row,1)+(N_seq_in_PAM+PAM_len+18-1);
                                variant_guide_exact_dsb_coord(row,1) = reference_sequence.coord(variant_guide_exact_dsb_index(row,1));
                                variant_guide_exact_strand(row,1) = '+';
                            elseif five_prime_PAM == 0 && variant_exact_match_list.strand(n,1) == 2 && (variant_exact_match_list.index(n,1)-PAM_len)>0
                                variant_guide_exact_coord_start(row,1) = reference_sequence.coord(variant_exact_match_list.index(n,1))-PAM_len;
                                variant_guide_exact_coord_end(row,1) = variant_guide_exact_coord_start(row,1)+sgRNA_len-1;
                                variant_guide_exact_index_start(row,1) = variant_exact_match_list.index(n,1)-PAM_len;
                                variant_guide_exact_index_end(row,1) = variant_exact_match_list.index(n,1)+sgRNA_len-1;
                                variant_guide_exact_dsb_index(row,1) = variant_guide_exact_index_start(row,1)+(PAM_len+2);
                                variant_guide_exact_dsb_coord(row,1) = reference_sequence.coord(variant_guide_exact_dsb_index(row,1));
                                variant_guide_exact_strand(row,1) = '-';
                            elseif five_prime_PAM == 1 && variant_exact_match_list.strand(n,1) == 2 && variant_exact_match_list.index(n,1) > 0
                                variant_guide_exact_coord_start(row,1) = reference_sequence.coord(variant_exact_match_list.index(n,1));
                                variant_guide_exact_coord_end(row,1) = reference_sequence.coord(variant_exact_match_list.index(n,1))+sgRNA_len+PAM_len-1;
                                variant_guide_exact_index_start(row,1) =variant_exact_match_list.index(n,1);
                                variant_guide_exact_index_end(row,1) = variant_exact_match_list.index(n,1)+PAM_len+sgRNA_len-1;
                                variant_guide_exact_dsb_index(row,1) = variant_guide_exact_index_start(row,1)+(N_seq_in_PAM+18);
                                variant_guide_exact_dsb_coord(row,1) = reference_sequence.coord(variant_guide_exact_dsb_index(row,1));
                                variant_guide_exact_strand(row,1) = '-';
                            end
                            row = row + 1;
                        end
                    end
                    
                    if length(find(variant_exact_match.index>0))>0;
                        for n = 1:length(variant_guide_exact_coord_start)
                            if strmatch(variant_guide_exact_strand(n,1),'+') == 1
                                if variant_guide_exact_index_end(n,1)>length(reference_seq);
                                    variant_guide_exact_index_end(n,1) = length(reference_seq);
                                end
                                variant_exact.guides{n,1} = strjoin(transpose(reference_sequence.nucleotide(variant_guide_exact_index_start(n,1):variant_guide_exact_index_end(n,1))),'');
                                variant_exact.strand{n,1} = '+';
                            elseif strmatch(variant_guide_exact_strand(n,1),'-') == 1
                                if variant_guide_exact_index_end(n,1)>length(reference_seq);
                                    variant_guide_exact_index_end(n,1) = length(reference_seq);
                                end
                                variant_exact.guides{n,1} = seqrcomplement(strjoin(transpose(reference_sequence.nucleotide(variant_guide_exact_index_start(n,1):variant_guide_exact_index_end(n,1))),''));
                                variant_exact.strand{n,1} = '-';
                            end
                            if five_prime_PAM == 0;
                                variant_exact.guides_only{n,1} = variant_exact.guides{n,1}(1:sgRNA_len);
                                variant_exact.pam_only{n,1} = variant_exact.guides{n,1}(sgRNA_len+1:length(variant_exact.guides{n,1}));
                            elseif five_prime_PAM == 1
                                variant_exact.guides_only{n,1} = variant_exact.guides{n,1}((PAM_len+1):(PAM_len+sgRNA_len));
                                variant_exact.pam_only{n,1} = variant_exact.guides{n,1}(1:PAM_len);
                            end
                            variant_exact.dsb_index(n,1) = variant_guide_exact_dsb_index(n,1);
                            variant_exact.dsb_coord(n,1) = variant_guide_exact_dsb_coord(n,1);
                            variant_exact.guide_index_start(n,1) = variant_guide_exact_index_start(n,1);
                            variant_exact.guide_index_end(n,1) = variant_guide_exact_index_end(n,1);
                            variant_exact.guide_coord_start(n,1) = variant_guide_exact_coord_start(n,1);
                            variant_exact.guide_coord_end(n,1) = variant_guide_exact_coord_end(n,1);
                        end
                    end
                    
                    unique_variant_guides_full.n_mm_wt = min_value;
                    if sum(min_index) > 0
                    unique_variant_guides_full.guides_wt = original_guides_full.guides(min_index);
                    unique_variant_guides_full.dsb_coord_wt = original_guides_full.dsb_coord(min_index);
                    unique_variant_guides_full.guides_only_wt = original_guides_full.guides_only(min_index);
                    unique_variant_guides_full.pam_only_wt = original_guides_full.pam_only(min_index);
                    unique_variant_guides_full.strand_wt = original_guides_full.strand(min_index);
                    unique_variant_guides_full.dsb_index_wt = original_guides_full.dsb_index(min_index);
                    unique_variant_guides_full.dsb_coord_wt = original_guides_full.dsb_coord(min_index);
                    unique_variant_guides_full.guide_index_start_wt = original_guides_full.guide_index_start(min_index);
                    unique_variant_guides_full.guide_index_end_wt = original_guides_full.guide_index_end(min_index);
                    unique_variant_guides_full.guide_coord_start_wt = original_guides_full.guide_coord_start(min_index);
                    unique_variant_guides_full.guide_coord_end_wt = original_guides_full.guide_coord_end(min_index);
                    end
                    if iscell(unique_variant_guides_full.n_mm_wt) == 0
                        unique_variant_guides_full.n_mm_wt = cellstr(num2str(unique_variant_guides_full.n_mm_wt));
                    end
                    
                    if length(find(variant_exact_match.index>0))>0
                        unique_variant_guides_full.n_mm_wt(find(variant_exact_match.index>0)) = {'N/A'};
                        unique_variant_guides_full.guides_wt(find(variant_exact_match.index>0)) = variant_exact.guides;
                        unique_variant_guides_full.dsb_coord_wt(find(variant_exact_match.index>0)) = variant_exact.dsb_coord;
                        unique_variant_guides_full.guides_only_wt(find(variant_exact_match.index>0)) = variant_exact.guides_only;
                        unique_variant_guides_full.pam_only_wt(find(variant_exact_match.index>0)) = variant_exact.pam_only;
                        unique_variant_guides_full.strand_wt(find(variant_exact_match.index>0)) = variant_exact.strand;
                        unique_variant_guides_full.dsb_index_wt(find(variant_exact_match.index>0)) = variant_exact.dsb_index;
                        unique_variant_guides_full.dsb_coord_wt(find(variant_exact_match.index>0)) = variant_exact.dsb_coord;
                        unique_variant_guides_full.guide_index_start_wt(find(variant_exact_match.index>0)) = variant_exact.guide_index_start;
                        unique_variant_guides_full.guide_index_end_wt(find(variant_exact_match.index>0)) = variant_exact.guide_index_end;
                        unique_variant_guides_full.guide_coord_start_wt(find(variant_exact_match.index>0)) = variant_exact.guide_coord_start;
                        unique_variant_guides_full.guide_coord_end_wt(find(variant_exact_match.index>0)) = variant_exact.guide_coord_end;
                        
                    end
                end
                
                chr_string = 'chr';
                if length(unique_variant_guides_full.type)>0 & length(original_guides_full.guides)>0
                    for m = 1:length(unique_variant_guides_full.type)
                        if strmatch(unique_variant_guides_full.strand(m,1),'+') == 1
                            idx = find((unique_variant_guides_full.guide_coord_end_wt(m)-test_variants.bpstart).*(test_variants.bpend-unique_variant_guides_full.guide_coord_start_wt(m))>=0);
                            num_mult_var_variant(m,1) = length(idx);
                            num_associated_variants(b,1) = max(num_mult_var_variant);
                            num_associated_variants_temp = num_associated_variants;
                            num_associated_variants_temp(b,:) = [];
                            num_spots(b,1) = length(num_mult_var_variant);
                            num_spots_sum(b,1) = sum(num_spots);
                            for e = 1:length(idx);
                                eval(['variant_associated_variants',num2str(e),'.bpstart(m,1) = cellstr(num2str(test_variants.bpstart(idx(e))));']);
                                eval(['variant_associated_variants',num2str(e),'.bpend(m,1) = cellstr(num2str(test_variants.bpend(idx(e))));']);
                                eval(['variant_associated_variants',num2str(e),'.chr_id(m,1) = test_variants.chr_id(idx(e));']);
                                eval(['variant_associated_variants',num2str(e),'.rs_id(m,1) = test_variants.rs_id(idx(e));']);
                                eval(['variant_associated_variants',num2str(e),'.ref_snp(m,1) = test_variants.ref_snp(idx(e));']);
                                eval(['variant_associated_variants',num2str(e),'.alt_snp(m,1) = test_variants.alt_snp(idx(e));']);
                                eval(['variant_associated_variants',num2str(e),'.identifier(m,1) = test_variants.type(idx(e));']);
                                eval(['variant_associated_variants',num2str(e),'.variant_index_start(m,1) = cellstr(num2str(test_variants.variant_index_start(idx(e))));']);
                                eval(['variant_associated_variants',num2str(e),'.variant_index_end(m,1) = cellstr(num2str(test_variants.variant_index_end(idx(e))));']);
                            end
                        elseif strmatch(unique_variant_guides_full.strand(m,1),'-') == 1
                            idx = find((unique_variant_guides_full.guide_coord_start_wt(m)-test_variants.bpstart).*(test_variants.bpend-unique_variant_guides_full.guide_coord_end_wt(m))>=0);
                            num_mult_var_variant(m,1) = length(idx);
                            num_associated_variants(b,1) = max(num_mult_var_variant);
                            num_associated_variants_temp = num_associated_variants;
                            num_associated_variants_temp(b,:) = [];
                            num_spots(b,1) = length(num_mult_var_variant);
                            num_spots_sum(b,1) = sum(num_spots);
                            for e = 1:length(idx);
                                eval(['variant_associated_variants',num2str(e),'.bpstart(m,1) = cellstr(num2str(test_variants.bpstart(idx(e))));']);
                                eval(['variant_associated_variants',num2str(e),'.bpend(m,1) = cellstr(num2str(test_variants.bpend(idx(e))));']);
                                eval(['variant_associated_variants',num2str(e),'.chr_id(m,1) = test_variants.chr_id(idx(e));']);
                                eval(['variant_associated_variants',num2str(e),'.rs_id(m,1) = test_variants.rs_id(idx(e));']);
                                eval(['variant_associated_variants',num2str(e),'.ref_snp(m,1) = test_variants.ref_snp(idx(e));']);
                                eval(['variant_associated_variants',num2str(e),'.alt_snp(m,1) = test_variants.alt_snp(idx(e));']);
                                eval(['variant_associated_variants',num2str(e),'.identifier(m,1) = test_variants.type(idx(e));']);
                                eval(['variant_associated_variants',num2str(e),'.variant_index_start(m,1) = cellstr(num2str(test_variants.variant_index_start(idx(e))));']);
                                eval(['variant_associated_variants',num2str(e),'.variant_index_end(m,1) = cellstr(num2str(test_variants.variant_index_end(idx(e))));']);
                            end
                        end
                    end
                    
                    for e = 1:max(num_mult_var_variant);
                        if eval(['length(variant_associated_variants',num2str(e),'.bpstart) < num_spots(b,1);'])
                            eval(['len_addition = length(variant_associated_variants',num2str(e),'.bpstart);']);
                            eval(['variant_associated_variants',num2str(e),'.bpstart((len_addition+1):num_spots(b,1),1) = {''''};'])
                            eval(['variant_associated_variants',num2str(e),'.bpend((len_addition+1):num_spots(b,1),1) = {''''};'])
                            eval(['variant_associated_variants',num2str(e),'.chr_id((len_addition+1):num_spots(b,1),1) = {''''};'])
                            eval(['variant_associated_variants',num2str(e),'.rs_id((len_addition+1):num_spots(b,1),1) = {''''};'])
                            eval(['variant_associated_variants',num2str(e),'.ref_snp((len_addition+1):num_spots(b,1),1) = {''''};'])
                            eval(['variant_associated_variants',num2str(e),'.alt_snp((len_addition+1):num_spots(b,1),1) = {''''};'])
                            eval(['variant_associated_variants',num2str(e),'.identifier((len_addition+1):num_spots(b,1),1) = {''''};'])
                            eval(['variant_associated_variants',num2str(e),'.variant_index_start((len_addition+1):num_spots(b,1),1) = {''''};'])
                            eval(['variant_associated_variants',num2str(e),'.variant_index_end((len_addition+1):num_spots(b,1),1) = {''''};'])
                        end
                    end
                end
                if twenty_extra_bases == 1 && perform_variant_analysis == 1;
                    edge_index_left = find(original_guides_full.dsb_coord<(seq_bpstart+sgRNA_len));
                    original_guides_full.guides(edge_index_left) = [];
                    original_guides_full.guides_only(edge_index_left) = [];
                    original_guides_full.pam_only(edge_index_left) = [];
                    original_guides_full.strand(edge_index_left) = [];
                    original_guides_full.dsb_index(edge_index_left) = [];
                    original_guides_full.dsb_coord(edge_index_left) = [];
                    original_guides_full.guide_index_start(edge_index_left) = [];
                    original_guides_full.guide_index_end(edge_index_left) = [];
                    original_guides_full.guide_coord_start(edge_index_left) = [];
                    original_guides_full.guide_coord_end(edge_index_left) = [];
                    edge_index_right = find(original_guides_full.dsb_coord>(seq_bpend-sgRNA_len));
                    original_guides_full.guides(edge_index_right) = [];
                    original_guides_full.guides_only(edge_index_right) = [];
                    original_guides_full.pam_only(edge_index_right) = [];
                    original_guides_full.strand(edge_index_right) = [];
                    original_guides_full.dsb_index(edge_index_right) = [];
                    original_guides_full.dsb_coord(edge_index_right) = [];
                    original_guides_full.guide_index_start(edge_index_right) = [];
                    original_guides_full.guide_index_end(edge_index_right) = [];
                    original_guides_full.guide_coord_start(edge_index_right) = [];
                    original_guides_full.guide_coord_end(edge_index_right) = [];
                end
            end
            if perform_variant_analysis == 1 && exist('original_guides_full','var')==1 || perform_wgs_analysis == 1 && exist('original_guides_full','var')==1;
                if a == starting_number;
                    for n = 1:length(original_guides_full.guides)
                        zeros_to_add = [];
                        if exist('guide_ID_num_temp','var') == 0
                            guide_ID_num_temp = n;
                        else guide_ID_num_temp = guide_ID_num_temp + 1;
                        end
                        num_zeros_add = length(num2str(guide_ID_num_temp));
                        if num_zeros_add < max_number_guide_ID_digits
                            for j = 1:(max_number_guide_ID_digits-num_zeros_add)
                                zeros_to_add = strcat(zeros_to_add,'0');
                            end
                            original_guides_full.guide_ID_number{n,1} = strcat('sgRNA_',zeros_to_add,num2str(guide_ID_num_temp));
                        elseif num_zeros_add == max_number_guide_ID_digits
                            original_guides_full.guide_ID_number{n,1} = ['sgRNA_',num2str(guide_ID_num_temp),''];
                        end
                    end
                elseif a > starting_number
                    for n = 1:length(original_guides_full.guides)
                        zeros_to_add = [];
                        guide_ID_num_temp = guide_ID_num_temp+1;
                        num_zeros_add = length(num2str(guide_ID_num_temp));
                        if num_zeros_add < max_number_guide_ID_digits
                            for j = 1:(max_number_guide_ID_digits-num_zeros_add)
                                zeros_to_add = strcat(zeros_to_add,'0');
                            end
                            original_guides_full.guide_ID_number{n,1} = strcat('sgRNA_',zeros_to_add,num2str(guide_ID_num_temp));
                        elseif num_zeros_add == max_number_guide_ID_digits
                            original_guides_full.guide_ID_number{n,1} = ['sgRNA_',num2str(guide_ID_num_temp),''];
                        end
                    end
                end
            end
          
            if twenty_extra_bases == 1 && perform_variant_analysis == 1 && exist('unique_variant_guides_full','var')==1 || twenty_extra_bases == 1 && perform_wgs_analysis == 1 && exist('unique_variant_guides_full','var')==1;
                if twenty_extra_bases == 1 && perform_variant_analysis == 1 && length(unique_variant_guides_full.guides)>0 && length(original_guides_full.guides)>0 || twenty_extra_bases == 1 && perform_wgs_analysis == 1 && length(unique_variant_guides_full.guides)>0 && length(original_guides_full.guides)>0;
                    edge_index_left = find(unique_variant_guides_full.dsb_coord_wt<(seq_bpstart+sgRNA_len));
                    unique_variant_guides_full.guides(edge_index_left) = [];
                    unique_variant_guides_full.guides_only(edge_index_left) = [];
                    unique_variant_guides_full.pam_only(edge_index_left) = [];
                    unique_variant_guides_full.strand(edge_index_left) = [];
                    unique_variant_guides_full.pam_create(edge_index_left) = [];
                    unique_variant_guides_full.n_mm_wt(edge_index_left) = [];
                    unique_variant_guides_full.guides_wt(edge_index_left) = [];
                    unique_variant_guides_full.guides_only_wt(edge_index_left) = [];
                    unique_variant_guides_full.pam_only_wt(edge_index_left) = [];
                    unique_variant_guides_full.strand_wt(edge_index_left) = [];
                    unique_variant_guides_full.dsb_index_wt(edge_index_left) = [];
                    unique_variant_guides_full.dsb_coord_wt(edge_index_left) = [];
                    unique_variant_guides_full.guide_index_start_wt(edge_index_left) = [];
                    unique_variant_guides_full.guide_index_end_wt(edge_index_left) = [];
                    unique_variant_guides_full.guide_coord_start_wt(edge_index_left) = [];
                    unique_variant_guides_full.guide_coord_end_wt(edge_index_left) = [];
                    
                    for e = 1:num_associated_variants(b,1)
                        eval(['variant_associated_variants',num2str(e),'.bpstart(edge_index_left) = [];']);
                        eval(['variant_associated_variants',num2str(e),'.bpend(edge_index_left) = [];']);
                        eval(['variant_associated_variants',num2str(e),'.rs_id(edge_index_left) = [];']);
                        eval(['variant_associated_variants',num2str(e),'.ref_snp(edge_index_left) = [];']);
                        eval(['variant_associated_variants',num2str(e),'.alt_snp(edge_index_left) = [];']);
                        eval(['variant_associated_variants',num2str(e),'.identifier(edge_index_left) = [];']);
                        eval(['variant_associated_variants',num2str(e),'.variant_index_start(edge_index_left) = [];']);
                        eval(['variant_associated_variants',num2str(e),'.variant_index_end(edge_index_left) = [];']);
                        eval(['variant_associated_variants',num2str(e),'.chr_id(edge_index_left) = [];']);
                    end
                    
                    edge_index_right = find(unique_variant_guides_full.dsb_coord_wt>(seq_bpend-sgRNA_len));
                    unique_variant_guides_full.guides(edge_index_right) = [];
                    unique_variant_guides_full.guides_only(edge_index_right) = [];
                    unique_variant_guides_full.pam_only(edge_index_right) = [];
                    unique_variant_guides_full.strand(edge_index_right) = [];
                    unique_variant_guides_full.pam_create(edge_index_right) = [];
                    
                    unique_variant_guides_full.n_mm_wt(edge_index_right) = [];
                    unique_variant_guides_full.guides_wt(edge_index_right) = [];
                    unique_variant_guides_full.guides_only_wt(edge_index_right) = [];
                    unique_variant_guides_full.pam_only_wt(edge_index_right) = [];
                    unique_variant_guides_full.strand_wt(edge_index_right) = [];
                    unique_variant_guides_full.dsb_index_wt(edge_index_right) = [];
                    unique_variant_guides_full.dsb_coord_wt(edge_index_right) = [];
                    unique_variant_guides_full.guide_index_start_wt(edge_index_right) = [];
                    unique_variant_guides_full.guide_index_end_wt(edge_index_right) = [];
                    unique_variant_guides_full.guide_coord_start_wt(edge_index_right) = [];
                    unique_variant_guides_full.guide_coord_end_wt(edge_index_right) = [];
                    
                    for e = 1:num_associated_variants(b,1)
                        eval(['variant_associated_variants',num2str(e),'.bpstart(edge_index_right) = [];']);
                        eval(['variant_associated_variants',num2str(e),'.bpend(edge_index_right) = [];']);
                        eval(['variant_associated_variants',num2str(e),'.rs_id(edge_index_right) = [];']);
                        eval(['variant_associated_variants',num2str(e),'.ref_snp(edge_index_right) = [];']);
                        eval(['variant_associated_variants',num2str(e),'.alt_snp(edge_index_right) = [];']);
                        eval(['variant_associated_variants',num2str(e),'.identifier(edge_index_right) = [];']);
                        eval(['variant_associated_variants',num2str(e),'.variant_index_start(edge_index_right) = [];']);
                        eval(['variant_associated_variants',num2str(e),'.variant_index_end(edge_index_right) = [];']);
                        eval(['variant_associated_variants',num2str(e),'.chr_id(edge_index_right) = [];']);
                    end
                    
                    num_spots(b,1) = num_spots(b,1)-length(edge_index_right)-length(edge_index_left);
                    num_spots_sum(b,1) = num_spots_sum(b,1)-length(edge_index_right)-length(edge_index_left);
                    
                end
            end
            
            if exist('num_associated_variants','var') == 1;
            num_associated_variants_batch(a,1) = max(num_associated_variants);
            end
            
            if exist('unique_variant_guides_full','var')==1
                if length(unique_variant_guides_full.guides)>0 && length(original_guides_full.guides)>0
                    for e = 1:length(unique_variant_guides_full.guides_only)
                        if length(strfind(reference_seq,char(unique_variant_guides_full.guides_only(e,1))))>0
                            unique_variant_guides_full.guide_ID_number_wt(e,1) = {'N/A'};
                        elseif length(strfind(reference_seq,seqrcomplement(char(unique_variant_guides_full.guides_only(e,1)))))>0
                            unique_variant_guides_full.guide_ID_number_wt(e,1) = {'N/A'};
                        elseif length(find(ismember(original_guides_full.guides_only,unique_variant_guides_full.guides_only_wt(e,1))>0))>0
                            unique_variant_guides_full.guide_ID_number_wt(e,1) = original_guides_full.guide_ID_number(find(ismember(original_guides_full.guides_only,unique_variant_guides_full.guides_only_wt(e,1))>0));
                        elseif length(find(ismember(original_guides_full.guides_only,unique_variant_guides_full.guides_only_wt(e,1))>0))>0
                            unique_variant_guides_full.guide_ID_number_wt(e,1) = original_guides_full.guide_ID_number(find(ismember(original_guides_full.guides_only,seqrcomplement(char(unique_variant_guides_full.guides_only_wt(e,1))))>0));
                        end
                    end
                    for n = 1:length(unique_variant_guides_full.guides)
                        zeros_to_add = [];
                        guide_ID_num_temp = guide_ID_num_temp+1;
                        num_zeros_add = length(num2str(guide_ID_num_temp));
                        if num_zeros_add < max_number_guide_ID_digits
                            for j = 1:(max_number_guide_ID_digits-num_zeros_add)
                                zeros_to_add = strcat(zeros_to_add,'0');
                            end
                            unique_variant_guides_full.guide_ID_number{n,1} = strcat('sgRNA_',zeros_to_add,num2str(guide_ID_num_temp));
                        elseif num_zeros_add == max_number_guide_ID_digits
                            unique_variant_guides_full.guide_ID_number{n,1} = ['sgRNA_',num2str(guide_ID_num_temp),''];
                        end
                    end
                end
            end
        end
        
        %% Output analysis
        
        if b == length(PAM_sequence_list)
            disp('Outputting analysis...')
            disp(' ')
        end
        
        temp_filename = ['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,''];
        output_filenames_list = [output_filenames_list;cellstr(temp_filename)];
        
        if exist('original_guides_full','var') == 1
            original_guides_full.chr_id(1:length(original_guides_full.guides),1) = ref_seq_info.chr(a,1);
            original_guides_full.pam(1:length(original_guides_full.guides),1) = PAM_sequence_list(b,1);
        end
        if multiple_match_analysis == 1 && exist('multiple_match','var')==1
            multiple_match.chr_id(1:length(multiple_match.guides),1) = ref_seq_info.chr(a,1);
            for u = 1:length(multiple_match.guides)
                if length(find(ismember(original_guides_full.guides_only,multiple_match.guides_only(u,1))>0)) > 1
                    for v = 1:length(find(ismember(original_guides_full.guides_only,multiple_match.guides_only(u,1))>0))
                        guide_num_temp_mult = original_guides_full.guide_ID_number(transpose(find(ismember(original_guides_full.guides_only,multiple_match.guides_only(u,1))>0)));
                        if v == length(find(ismember(original_guides_full.guides_only,multiple_match.guides_only(u,1))>0))
                            multiple_match.guide_ID_number(u,1) = cellstr(strjoin(guide_num_temp_mult));
                        end
                    end
                else
                    multiple_match.guide_ID_number(u,1) = original_guides_full.guide_ID_number(find(ismember(original_guides_full.guides_only,multiple_match.guides_only(u,1))>0));
                end
            end
        end
        if perform_variant_analysis == 1 && exist('unique_variant_guides_full','var')==1 || perform_wgs_analysis == 1 && exist('unique_variant_guides_full','var')==1
            unique_variant_guides_full.chr_id_wt(1:length(unique_variant_guides_full.guides),1) = ref_seq_info.chr(a,1);
            unique_variant_guides_full.pam(1:length(unique_variant_guides_full.guides),1) = {PAM_sequence1};
            unique_variant_guides_full.chr_id(1:length(unique_variant_guides_full.guides),1) = ref_seq_info.chr(a,1);
        end
        if perform_haplotype_analysis == 1
            unique_haplotype_guides_full.chr_id_wt(1:length(unique_haplotype_guides_full.guides),1) = ref_seq_info.chr(a,1);
            unique_haplotype_guides_full.chr_id(1:length(unique_haplotype_guides_full.guides),1) = ref_seq_info.chr(a,1);
            unique_haplotype_guides_full.pam(1:length(unique_haplotype_guides_full.guides),1) = PAM_sequence_list(b,1);
        end
        temp_coord_string_start = cellstr(strcat('Guide Start (',coordinate_sys,' Coordinate)'));
        temp_coord_string_end = cellstr(strcat('Guide End (',coordinate_sys,' Coordinate)'));
        temp_coord_string_start_txt = cellstr(strcat('Guide_Start_',coordinate_sys,'_Coordinate'));
        temp_coord_string_end_txt = cellstr(strcat('Guide_End_',coordinate_sys,'_Coordinate'));
        
        temp_coord_string_start_variant = cellstr(strcat('Variant Start (',coordinate_sys,' Coordinate)'));
        temp_coord_string_end_variant = cellstr(strcat('Variant End (',coordinate_sys,' Coordinate)'));
        temp_coord_string_dsb = cellstr(strcat('Double Strand Break (',coordinate_sys,' Coordinate)'));
        temp_coord_string_dsb_txt = cellstr(strcat('Double_Strand_Break_',coordinate_sys,'_Coordinate'));
        
        if length(cellstr(PAM_sequence1)) > 1
            PAM_sequence1 = strjoin(PAM_sequence1);
        end
        
        if exist('original_guides_full','var')==1
        if length(original_guides_full.guides)>0
            merge_output_filenames_list_temp(1:length(original_guides_full.guides)) = cellstr(temp_filename);
            merge_original_guides_full.guide_ID_number = [merge_original_guides_full.guide_ID_number;original_guides_full.guide_ID_number];
            merge_original_guides_full.pam = [merge_original_guides_full.pam;original_guides_full.pam];
            merge_original_guides_full.guides_only = [merge_original_guides_full.guides_only;original_guides_full.guides_only];
            merge_original_guides_full.pam_only = [merge_original_guides_full.pam_only;original_guides_full.pam_only];
            merge_original_guides_full.guides = [merge_original_guides_full.guides;original_guides_full.guides];
            merge_original_guides_full.chr_id = [merge_original_guides_full.chr_id;original_guides_full.chr_id];
            merge_original_guides_full.strand = [merge_original_guides_full.strand;original_guides_full.strand];
            merge_original_guides_full.guide_coord_start = [merge_original_guides_full.guide_coord_start;original_guides_full.guide_coord_start];
            merge_original_guides_full.guide_coord_end = [merge_original_guides_full.guide_coord_end;original_guides_full.guide_coord_end];
            merge_original_guides_full.dsb_coord = [merge_original_guides_full.dsb_coord;original_guides_full.dsb_coord];
            merge_original_guides_full.guide_index_start = [merge_original_guides_full.guide_index_start;original_guides_full.guide_index_start];
            merge_original_guides_full.guide_index_end = [merge_original_guides_full.guide_index_end;original_guides_full.guide_index_end];
            merge_original_guides_full.dsb_index = [merge_original_guides_full.dsb_index;original_guides_full.dsb_index];
            if multiple_match_analysis == 1
                merge_original_guides_full.mult_match_count = [merge_original_guides_full.mult_match_count;original_guides_full.mult_match_count];
            end
            merge_original_guides_full.filename = [merge_original_guides_full.filename;transpose(merge_output_filenames_list_temp)];
        end
        end
        if end_number > starting_number && perform_batch_output == 1 && b == length(PAM_sequence_list) && exist('merge_original_guides_full','var')==1
            batch_output_filenames_list_temp(1:length(merge_original_guides_full.guides)) = cellstr(temp_filename);
            batch_original_guides_full.guide_ID_number = [batch_original_guides_full.guide_ID_number;merge_original_guides_full.guide_ID_number];
            batch_original_guides_full.pam = [batch_original_guides_full.pam;merge_original_guides_full.pam];
            batch_original_guides_full.guides_only = [batch_original_guides_full.guides_only;merge_original_guides_full.guides_only];
            batch_original_guides_full.pam_only = [batch_original_guides_full.pam_only;merge_original_guides_full.pam_only];
            batch_original_guides_full.guides = [batch_original_guides_full.guides;merge_original_guides_full.guides];
            batch_original_guides_full.chr_id = [batch_original_guides_full.chr_id;merge_original_guides_full.chr_id];
            batch_original_guides_full.strand = [batch_original_guides_full.strand;merge_original_guides_full.strand];
            batch_original_guides_full.guide_coord_start = [batch_original_guides_full.guide_coord_start;merge_original_guides_full.guide_coord_start];
            batch_original_guides_full.guide_coord_end = [batch_original_guides_full.guide_coord_end;merge_original_guides_full.guide_coord_end];
            batch_original_guides_full.dsb_coord = [batch_original_guides_full.dsb_coord;merge_original_guides_full.dsb_coord];
            batch_original_guides_full.guide_index_start = [batch_original_guides_full.guide_index_start;merge_original_guides_full.guide_index_start];
            batch_original_guides_full.guide_index_end = [batch_original_guides_full.guide_index_end;merge_original_guides_full.guide_index_end];
            batch_original_guides_full.dsb_index = [batch_original_guides_full.dsb_index;merge_original_guides_full.dsb_index];
            batch_original_guides_full.mult_match_count = [batch_original_guides_full.mult_match_count;merge_original_guides_full.mult_match_count];
            batch_original_guides_full.filename = [batch_original_guides_full.filename;transpose(batch_output_filenames_list_temp)];
        end
        
        if perform_batch_output == 1 && exist('multiple_match','var') == 1
            merge_multiple_match.guide_ID_number = [merge_multiple_match.guide_ID_number;multiple_match.guide_ID_number];
            merge_multiple_match.pam = [merge_multiple_match.pam;multiple_match.pam];
            merge_multiple_match.guides_only = [merge_multiple_match.guides_only;multiple_match.guides_only];
            merge_multiple_match.pam_only = [merge_multiple_match.pam_only;multiple_match.pam_only];
            merge_multiple_match.guides = [merge_multiple_match.guides;multiple_match.guides];
            merge_multiple_match.strand = [merge_multiple_match.strand;multiple_match.strand];
            merge_multiple_match.chr_id = [merge_multiple_match.chr_id;multiple_match.chr_id];
            merge_multiple_match.guide_coord_start = [merge_multiple_match.guide_coord_start;multiple_match.guide_coord_start];
            merge_multiple_match.guide_coord_end = [merge_multiple_match.guide_coord_end;multiple_match.guide_coord_end];
            merge_multiple_match.dsb_coord = [merge_multiple_match.dsb_coord;multiple_match.dsb_coord];
            merge_multiple_match.guide_index_start = [merge_multiple_match.guide_index_start;multiple_match.guide_index_start];
            merge_multiple_match.guide_index_end = [merge_multiple_match.guide_index_end;multiple_match.guide_index_end];
            merge_multiple_match.dsb_index = [merge_multiple_match.dsb_index;multiple_match.dsb_index];
            merge_output_filenames_list_temp_mult(1:length(multiple_match.guides)) = cellstr(temp_filename);
            merge_multiple_match.filename = [merge_multiple_match.filename;transpose(merge_output_filenames_list_temp_mult)];
        end
        if end_number > starting_number && perform_batch_output == 1 && multiple_match_analysis == 1 && exist('merge_multiple_match','var') == 1 && b == length(PAM_sequence_list)
            batch_multiple_match.guide_ID_number = [batch_multiple_match.guide_ID_number;merge_multiple_match.guide_ID_number];
            batch_multiple_match.pam = [batch_multiple_match.pam;merge_multiple_match.pam];
            batch_multiple_match.guides_only = [batch_multiple_match.guides_only;merge_multiple_match.guides_only];
            batch_multiple_match.pam_only = [batch_multiple_match.pam_only;merge_multiple_match.pam_only];
            batch_multiple_match.guides = [batch_multiple_match.guides;merge_multiple_match.guides];
            batch_multiple_match.strand = [batch_multiple_match.strand;merge_multiple_match.strand];
            batch_multiple_match.chr_id = [batch_multiple_match.chr_id;merge_multiple_match.chr_id];
            batch_multiple_match.guide_coord_start = [batch_multiple_match.guide_coord_start;merge_multiple_match.guide_coord_start];
            batch_multiple_match.guide_coord_end = [batch_multiple_match.guide_coord_end;merge_multiple_match.guide_coord_end];
            batch_multiple_match.dsb_coord = [batch_multiple_match.dsb_coord;merge_multiple_match.dsb_coord];
            batch_multiple_match.guide_index_start = [batch_multiple_match.guide_index_start;merge_multiple_match.guide_index_start];
            batch_multiple_match.guide_index_end = [batch_multiple_match.guide_index_end;merge_multiple_match.guide_index_end];
            batch_multiple_match.dsb_index = [batch_multiple_match.dsb_index;merge_multiple_match.dsb_index];
            batch_output_filenames_list_temp_mult(1:length(merge_multiple_match.guides)) = cellstr(temp_filename);
            batch_multiple_match.filename = [batch_multiple_match.filename;transpose(batch_output_filenames_list_temp_mult)];
        end
               
        if perform_variant_analysis == 1 && exist('unique_variant_guides_full','var')==1 || perform_wgs_analysis == 1 && exist('unique_variant_guides_full','var')==1
            if length(unique_variant_guides_full.guides)>0
                merge_output_filenames_list_variant_temp(1:length(unique_variant_guides_full.guides)) = cellstr(temp_filename);
                if exist('unique_variant_guides_full.guide_ID_number','var') == 1
                    merge_unique_variant_guides_full.guide_ID_number = [merge_unique_variant_guides_full.guide_ID_number;unique_variant_guides_full.guide_ID_number];
                elseif exist('unique_variant_guides_full.guide_ID_number','var') == 0
                    unique_variant_guides_full.guide_ID_number = [];
                    merge_unique_variant_guides_full.guide_ID_number = [merge_unique_variant_guides_full.guide_ID_number;unique_variant_guides_full.guide_ID_number];
                end
                merge_unique_variant_guides_full.guides_only = [merge_unique_variant_guides_full.guides_only;unique_variant_guides_full.guides_only];
                merge_unique_variant_guides_full.pam = [merge_unique_variant_guides_full.pam;unique_variant_guides_full.pam];
                merge_unique_variant_guides_full.pam_only = [merge_unique_variant_guides_full.pam_only;unique_variant_guides_full.pam_only];
                merge_unique_variant_guides_full.guides = [merge_unique_variant_guides_full.guides;unique_variant_guides_full.guides];
                merge_unique_variant_guides_full.chr_id = [merge_unique_variant_guides_full.chr_id;unique_variant_guides_full.chr_id];
                merge_unique_variant_guides_full.strand = [merge_unique_variant_guides_full.strand;unique_variant_guides_full.strand];
                merge_unique_variant_guides_full.filename = [merge_unique_variant_guides_full.filename;transpose(merge_output_filenames_list_variant_temp)];                
                merge_unique_variant_guides_full.pam_create = [merge_unique_variant_guides_full.pam_create;unique_variant_guides_full.pam_create];
                if length(original_guides_full.guides)>0
                merge_unique_variant_guides_full.n_mm_wt = [merge_unique_variant_guides_full.n_mm_wt;unique_variant_guides_full.n_mm_wt];
                merge_unique_variant_guides_full.guides_wt = [merge_unique_variant_guides_full.guides_wt;unique_variant_guides_full.guides_wt];
                merge_unique_variant_guides_full.guides_only_wt = [merge_unique_variant_guides_full.guides_only_wt;unique_variant_guides_full.guides_only_wt];
                merge_unique_variant_guides_full.guide_ID_number_wt = [merge_unique_variant_guides_full.guide_ID_number_wt;unique_variant_guides_full.guide_ID_number_wt];
                merge_unique_variant_guides_full.chr_id_wt = [merge_unique_variant_guides_full.chr_id_wt;unique_variant_guides_full.chr_id_wt];
                merge_unique_variant_guides_full.pam_only_wt = [merge_unique_variant_guides_full.pam_only_wt;unique_variant_guides_full.pam_only_wt];
                merge_unique_variant_guides_full.strand_wt = [merge_unique_variant_guides_full.strand_wt;unique_variant_guides_full.strand_wt];
                merge_unique_variant_guides_full.dsb_index_wt = [merge_unique_variant_guides_full.dsb_index_wt;unique_variant_guides_full.dsb_index_wt];
                merge_unique_variant_guides_full.dsb_coord_wt = [merge_unique_variant_guides_full.dsb_coord_wt;unique_variant_guides_full.dsb_coord_wt];
                merge_unique_variant_guides_full.guide_index_start_wt = [merge_unique_variant_guides_full.guide_index_start_wt;unique_variant_guides_full.guide_index_start_wt];
                merge_unique_variant_guides_full.guide_index_end_wt = [merge_unique_variant_guides_full.guide_index_end_wt;unique_variant_guides_full.guide_index_end_wt];
                merge_unique_variant_guides_full.guide_coord_start_wt = [merge_unique_variant_guides_full.guide_coord_start_wt;unique_variant_guides_full.guide_coord_start_wt];
                merge_unique_variant_guides_full.guide_coord_end_wt = [merge_unique_variant_guides_full.guide_coord_end_wt;unique_variant_guides_full.guide_coord_end_wt];
                
                num_iterate = num_associated_variants(b,1);
                if max(num_associated_variants_temp) > num_iterate
                    num_iterate = max(num_associated_variants_temp);
                end
                for h = 1:num_iterate;
                    if eval(['exist([''merge_unique_variant_guides_full',num2str(h),'''],''var'')==0'])
                        
                        eval(['merge_unique_variant_guides_full',num2str(h),'.bpstart_var = [];']);
                        eval(['merge_unique_variant_guides_full',num2str(h),'.bpend_var = [];']);
                        eval(['merge_unique_variant_guides_full',num2str(h),'.rs_id_var = [];']);
                        eval(['merge_unique_variant_guides_full',num2str(h),'.ref_snp_var = [];']);
                        eval(['merge_unique_variant_guides_full',num2str(h),'.alt_snp_var = [];']);
                        eval(['merge_unique_variant_guides_full',num2str(h),'.identifier_var = [];']);
                        eval(['merge_unique_variant_guides_full',num2str(h),'.variant_index_start_var = [];']);
                        eval(['merge_unique_variant_guides_full',num2str(h),'.variant_index_end_var = [];']);
                        eval(['merge_unique_variant_guides_full',num2str(h),'.chr_id_var = [];']);
                    end
                end
                
                num_iterate = num_associated_variants(b,1);
                if max(num_associated_variants_temp) > num_iterate
                    num_iterate = max(num_associated_variants_temp);
                    for j = (num_associated_variants(b,1)+1):max(num_associated_variants_temp)
                        if eval(['exist([''variant_associated_variants',num2str(j),'''],''var'')==0'])
                            eval(['variant_associated_variants',num2str(j),'.bpstart(1:num_spots(b,1),1) = {''''};']);
                            eval(['variant_associated_variants',num2str(j),'.bpend(1:num_spots(b,1),1) = {''''};']);
                            eval(['variant_associated_variants',num2str(j),'.rs_id(1:num_spots(b,1),1) = {''''};']);
                            eval(['variant_associated_variants',num2str(j),'.ref_snp(1:num_spots(b,1),1) = {''''};']);
                            eval(['variant_associated_variants',num2str(j),'.alt_snp(1:num_spots(b,1),1) = {''''};']);
                            eval(['variant_associated_variants',num2str(j),'.identifier(1:num_spots(b,1),1) = {''''};']);
                            eval(['variant_associated_variants',num2str(j),'.variant_index_start(1:num_spots(b,1),1) = {''''};']);
                            eval(['variant_associated_variants',num2str(j),'.variant_index_end(1:num_spots(b,1),1) = {''''};']);
                            eval(['variant_associated_variants',num2str(j),'.chr_id(1:num_spots(b,1),1) = {''''};']);
                        end
                    end
                end
                if num_iterate > max(num_associated_variants_temp)
                    for j = (max(num_associated_variants_temp)+1):num_associated_variants(b,1)
                        eval(['variant_associated_variants',num2str(j),'.bpstart((num_spots_sum(b-1)+1):num_spots_sum(b,1),1) = variant_associated_variants',num2str(h),'.bpstart;']);
                        eval(['variant_associated_variants',num2str(j),'.bpend((num_spots_sum(b-1)+1):num_spots_sum(b,1),1) = variant_associated_variants',num2str(h),'.bpend;']);
                        eval(['variant_associated_variants',num2str(j),'.rs_id((num_spots_sum(b-1)+1):num_spots_sum(b,1),1) = variant_associated_variants',num2str(h),'.rs_id;']);
                        eval(['variant_associated_variants',num2str(j),'.ref_snp((num_spots_sum(b-1)+1):num_spots_sum(b,1),1) = variant_associated_variants',num2str(h),'.ref_snp;']);
                        eval(['variant_associated_variants',num2str(j),'.alt_snp((num_spots_sum(b-1)+1):num_spots_sum(b,1),1) = variant_associated_variants',num2str(h),'.alt_snp;']);
                        eval(['variant_associated_variants',num2str(j),'.identifier((num_spots_sum(b-1)+1):num_spots_sum(b,1),1) = variant_associated_variants',num2str(h),'.identifier;']);
                        eval(['variant_associated_variants',num2str(j),'.variant_index_start((num_spots_sum(b-1)+1):num_spots_sum(b,1),1) = variant_associated_variants',num2str(h),'.variant_index_start;']);
                        eval(['variant_associated_variants',num2str(j),'.variant_index_end((num_spots_sum(b-1)+1):num_spots_sum(b,1),1) = variant_associated_variants',num2str(h),'.variant_index_end;']);
                        eval(['variant_associated_variants',num2str(j),'.chr_id((num_spots_sum(b-1)+1):num_spots_sum(b,1),1) = variant_associated_variants',num2str(h),'.chr_id;']);
                        
                        eval(['variant_associated_variants',num2str(j),'.bpstart(1:num_spots(b,1),1) = {''''};']);
                        eval(['variant_associated_variants',num2str(j),'.bpend(1:num_spots(b,1),1) = {''''};']);
                        eval(['variant_associated_variants',num2str(j),'.rs_id(1:num_spots(b,1),1) = {''''};']);
                        eval(['variant_associated_variants',num2str(j),'.ref_snp(1:num_spots(b,1),1) = {''''};']);
                        eval(['variant_associated_variants',num2str(j),'.alt_snp(1:num_spots(b,1),1) = {''''};']);
                        eval(['variant_associated_variants',num2str(j),'.identifier(1:num_spots(b,1),1) = {''''};']);
                        eval(['variant_associated_variants',num2str(j),'.variant_index_start(1:num_spots(b,1),1) = {''''};']);
                        eval(['variant_associated_variants',num2str(j),'.variant_index_end(1:num_spots(b,1),1) = {''''};']);
                        eval(['variant_associated_variants',num2str(j),'.chr_id(1:num_spots(b,1),1) = {''''};']);
                    end
                end
                
                for h = 1:num_iterate;
                    eval(['merge_unique_variant_guides_full',num2str(h),'.bpstart_var = [merge_unique_variant_guides_full',num2str(h),'.bpstart_var;variant_associated_variants',num2str(h),'.bpstart];']);
                    eval(['merge_unique_variant_guides_full',num2str(h),'.bpend_var = [merge_unique_variant_guides_full',num2str(h),'.bpend_var;variant_associated_variants',num2str(h),'.bpend];']);
                    eval(['merge_unique_variant_guides_full',num2str(h),'.rs_id_var = [merge_unique_variant_guides_full',num2str(h),'.rs_id_var;variant_associated_variants',num2str(h),'.rs_id];']);
                    eval(['merge_unique_variant_guides_full',num2str(h),'.ref_snp_var = [merge_unique_variant_guides_full',num2str(h),'.ref_snp_var;variant_associated_variants',num2str(h),'.ref_snp];']);
                    eval(['merge_unique_variant_guides_full',num2str(h),'.alt_snp_var = [merge_unique_variant_guides_full',num2str(h),'.alt_snp_var;variant_associated_variants',num2str(h),'.alt_snp];']);
                    eval(['merge_unique_variant_guides_full',num2str(h),'.identifier_var = [merge_unique_variant_guides_full',num2str(h),'.identifier_var;variant_associated_variants',num2str(h),'.identifier];']);
                    eval(['merge_unique_variant_guides_full',num2str(h),'.variant_index_start_var = [merge_unique_variant_guides_full',num2str(h),'.variant_index_start_var;variant_associated_variants',num2str(h),'.variant_index_start];']);
                    eval(['merge_unique_variant_guides_full',num2str(h),'.variant_index_end_var = [merge_unique_variant_guides_full',num2str(h),'.variant_index_end_var;variant_associated_variants',num2str(h),'.variant_index_end];']);
                    eval(['merge_unique_variant_guides_full',num2str(h),'.chr_id_var = [merge_unique_variant_guides_full',num2str(h),'.chr_id_var;variant_associated_variants',num2str(h),'.chr_id];']);
                end
                
                if b == length(PAM_sequence_list)
                    num_spots_batch(a,1) = sum(num_spots);
                    num_spots_batch_temp = num_spots_batch;
                    num_spots_batch_temp(a,:) = [];
                    if length(num_spots_batch) == 1
                        num_spots_batch_temp = 0;
                    elseif length(num_spots_batch)>1
                        num_spots_batch_temp = sum(num_spots_batch_temp);
                    end
                end
                end
            end
        end
        if end_number > starting_number && perform_batch_output == 1 && perform_variant_analysis == 1 && length(merge_unique_variant_guides_full.guide_index_start_wt)>0 && b == length(PAM_sequence_list) || end_number > starting_number && perform_batch_output == 1 && perform_wgs_analysis == 1 && length(merge_unique_variant_guides_full.guide_index_start_wt)>0 && b == length(PAM_sequence_list)
            batch_output_filenames_list_variant_temp(1:length(merge_unique_variant_guides_full.guides)) = cellstr(temp_filename);
            batch_unique_variant_guides_full.guide_ID_number = [batch_unique_variant_guides_full.guide_ID_number;merge_unique_variant_guides_full.guide_ID_number];
            batch_unique_variant_guides_full.guides_only = [batch_unique_variant_guides_full.guides_only;merge_unique_variant_guides_full.guides_only];
            batch_unique_variant_guides_full.pam = [batch_unique_variant_guides_full.pam;merge_unique_variant_guides_full.pam];
            batch_unique_variant_guides_full.pam_only = [batch_unique_variant_guides_full.pam_only;merge_unique_variant_guides_full.pam_only];
            batch_unique_variant_guides_full.guides = [batch_unique_variant_guides_full.guides;merge_unique_variant_guides_full.guides];
            batch_unique_variant_guides_full.chr_id = [batch_unique_variant_guides_full.chr_id;merge_unique_variant_guides_full.chr_id];
            batch_unique_variant_guides_full.strand = [batch_unique_variant_guides_full.strand;merge_unique_variant_guides_full.strand];
            batch_unique_variant_guides_full.filename = [batch_unique_variant_guides_full.filename;transpose(batch_output_filenames_list_variant_temp)];
            
            batch_unique_variant_guides_full.pam_create = [batch_unique_variant_guides_full.pam_create;merge_unique_variant_guides_full.pam_create];
            batch_unique_variant_guides_full.n_mm_wt = [batch_unique_variant_guides_full.n_mm_wt;merge_unique_variant_guides_full.n_mm_wt];
            batch_unique_variant_guides_full.guides_wt = [batch_unique_variant_guides_full.guides_wt;merge_unique_variant_guides_full.guides_wt];
            batch_unique_variant_guides_full.guides_only_wt = [batch_unique_variant_guides_full.guides_only_wt;merge_unique_variant_guides_full.guides_only_wt];
            batch_unique_variant_guides_full.guide_ID_number_wt = [batch_unique_variant_guides_full.guide_ID_number_wt;merge_unique_variant_guides_full.guide_ID_number_wt];
            batch_unique_variant_guides_full.chr_id_wt = [batch_unique_variant_guides_full.chr_id_wt;merge_unique_variant_guides_full.chr_id_wt];
            batch_unique_variant_guides_full.pam_only_wt = [batch_unique_variant_guides_full.pam_only_wt;merge_unique_variant_guides_full.pam_only_wt];
            batch_unique_variant_guides_full.strand_wt = [batch_unique_variant_guides_full.strand_wt;merge_unique_variant_guides_full.strand_wt];
            batch_unique_variant_guides_full.dsb_index_wt = [batch_unique_variant_guides_full.dsb_index_wt;merge_unique_variant_guides_full.dsb_index_wt];
            batch_unique_variant_guides_full.dsb_coord_wt = [batch_unique_variant_guides_full.dsb_coord_wt;merge_unique_variant_guides_full.dsb_coord_wt];
            batch_unique_variant_guides_full.guide_index_start_wt = [batch_unique_variant_guides_full.guide_index_start_wt;merge_unique_variant_guides_full.guide_index_start_wt];
            batch_unique_variant_guides_full.guide_index_end_wt = [batch_unique_variant_guides_full.guide_index_end_wt;merge_unique_variant_guides_full.guide_index_end_wt];
            batch_unique_variant_guides_full.guide_coord_start_wt = [batch_unique_variant_guides_full.guide_coord_start_wt;merge_unique_variant_guides_full.guide_coord_start_wt];
            batch_unique_variant_guides_full.guide_coord_end_wt = [batch_unique_variant_guides_full.guide_coord_end_wt;merge_unique_variant_guides_full.guide_coord_end_wt];
            
            num_iterate = num_associated_variants_batch(a,1);
            num_associated_variants_batch_temp = num_associated_variants_batch;
            if a > 1
                num_associated_variants_batch_temp(a,:) = [];
            end
            
            for h = 1:num_iterate;
                if eval(['exist([''batch_unique_variant_guides_full',num2str(h),'''],''var'')==0'])
                    eval(['batch_unique_variant_guides_full',num2str(h),'.bpstart_var = {};']);
                    eval(['batch_unique_variant_guides_full',num2str(h),'.bpend_var = {};']);
                    eval(['batch_unique_variant_guides_full',num2str(h),'.rs_id_var = {};']);
                    eval(['batch_unique_variant_guides_full',num2str(h),'.ref_snp_var = {};']);
                    eval(['batch_unique_variant_guides_full',num2str(h),'.alt_snp_var = {};']);
                    eval(['batch_unique_variant_guides_full',num2str(h),'.identifier_var = {};']);
                    eval(['batch_unique_variant_guides_full',num2str(h),'.variant_index_start_var = {};']);
                    eval(['batch_unique_variant_guides_full',num2str(h),'.variant_index_end_var = {};']);
                    eval(['batch_unique_variant_guides_full',num2str(h),'.chr_id_var = {};']);
                end
            end
            
            num_iterate = num_associated_variants_batch(a,1);
            if max(num_associated_variants_batch_temp) > num_iterate
                for j = (num_associated_variants_batch(a,1)+1):max(num_associated_variants_batch_temp)
                    eval(['merge_unique_variant_guides_full',num2str(j),'.bpstart_var(1:num_spots_batch(a,1),1) = {''''};']);
                    eval(['merge_unique_variant_guides_full',num2str(j),'.bpend_var(1:num_spots_batch(a,1),1) = {''''};']);
                    eval(['merge_unique_variant_guides_full',num2str(j),'.rs_id_var(1:num_spots_batch(a,1),1) = {''''};']);
                    eval(['merge_unique_variant_guides_full',num2str(j),'.ref_snp_var(1:num_spots_batch(a,1),1) = {''''};']);
                    eval(['merge_unique_variant_guides_full',num2str(j),'.alt_snp_var(1:num_spots_batch(a,1),1) = {''''};']);
                    eval(['merge_unique_variant_guides_full',num2str(j),'.identifier_var(1:num_spots_batch(a,1),1) = {''''};']);
                    eval(['merge_unique_variant_guides_full',num2str(j),'.variant_index_start_var(1:num_spots_batch(a,1),1) = {''''};']);
                    eval(['merge_unique_variant_guides_full',num2str(j),'.variant_index_end_var(1:num_spots_batch(a,1),1) = {''''};']);
                    eval(['merge_unique_variant_guides_full',num2str(j),'.chr_id_var(1:num_spots_batch(a,1),1) = {''''};']);
                end
                for h = 1:max(num_associated_variants_batch_temp);
                    eval(['batch_unique_variant_guides_full',num2str(h),'.bpstart_var = [batch_unique_variant_guides_full',num2str(h),'.bpstart_var;merge_unique_variant_guides_full',num2str(h),'.bpstart_var];']);
                    eval(['batch_unique_variant_guides_full',num2str(h),'.bpend_var = [batch_unique_variant_guides_full',num2str(h),'.bpend_var;merge_unique_variant_guides_full',num2str(h),'.bpend_var];']);
                    eval(['batch_unique_variant_guides_full',num2str(h),'.rs_id_var = [batch_unique_variant_guides_full',num2str(h),'.rs_id_var;merge_unique_variant_guides_full',num2str(h),'.rs_id_var];']);
                    eval(['batch_unique_variant_guides_full',num2str(h),'.ref_snp_var = [batch_unique_variant_guides_full',num2str(h),'.ref_snp_var;merge_unique_variant_guides_full',num2str(h),'.ref_snp_var];']);
                    eval(['batch_unique_variant_guides_full',num2str(h),'.alt_snp_var = [batch_unique_variant_guides_full',num2str(h),'.alt_snp_var;merge_unique_variant_guides_full',num2str(h),'.alt_snp_var];']);
                    eval(['batch_unique_variant_guides_full',num2str(h),'.identifier_var = [batch_unique_variant_guides_full',num2str(h),'.identifier_var;merge_unique_variant_guides_full',num2str(h),'.identifier_var];']);
                    eval(['batch_unique_variant_guides_full',num2str(h),'.variant_index_start_var = [batch_unique_variant_guides_full',num2str(h),'.variant_index_start_var;merge_unique_variant_guides_full',num2str(h),'.variant_index_start_var];']);
                    eval(['batch_unique_variant_guides_full',num2str(h),'.variant_index_end_var = [batch_unique_variant_guides_full',num2str(h),'.variant_index_end_var;merge_unique_variant_guides_full',num2str(h),'.variant_index_end_var];']);
                    eval(['batch_unique_variant_guides_full',num2str(h),'.chr_id_var = [batch_unique_variant_guides_full',num2str(h),'.chr_id_var;merge_unique_variant_guides_full',num2str(h),'.chr_id_var];']);
                end
            elseif num_iterate > max(num_associated_variants_batch_temp)
                for j = 1:max(num_associated_variants_batch(a,1))
                    eval(['batch_unique_variant_guides_full',num2str(j),'.bpstart_var((num_spots_batch_temp+1):sum(num_spots_batch),1) = merge_unique_variant_guides_full',num2str(j),'.bpstart_var;']);
                    eval(['batch_unique_variant_guides_full',num2str(j),'.bpend_var((num_spots_batch_temp+1):sum(num_spots_batch),1) = merge_unique_variant_guides_full',num2str(j),'.bpend_var;']);
                    eval(['batch_unique_variant_guides_full',num2str(j),'.rs_id_var((num_spots_batch_temp+1):sum(num_spots_batch),1) = merge_unique_variant_guides_full',num2str(j),'.rs_id_var;']);
                    eval(['batch_unique_variant_guides_full',num2str(j),'.ref_snp_var((num_spots_batch_temp+1):sum(num_spots_batch),1) = merge_unique_variant_guides_full',num2str(j),'.ref_snp_var;']);
                    eval(['batch_unique_variant_guides_full',num2str(j),'.alt_snp_var((num_spots_batch_temp+1):sum(num_spots_batch),1) = merge_unique_variant_guides_full',num2str(j),'.alt_snp_var;']);
                    eval(['batch_unique_variant_guides_full',num2str(j),'.identifier_var((num_spots_batch_temp+1):sum(num_spots_batch),1) = merge_unique_variant_guides_full',num2str(j),'.identifier_var;']);
                    eval(['batch_unique_variant_guides_full',num2str(j),'.variant_index_start_var((num_spots_batch_temp+1):sum(num_spots_batch),1) = merge_unique_variant_guides_full',num2str(j),'.variant_index_start_var;']);
                    eval(['batch_unique_variant_guides_full',num2str(j),'.variant_index_end_var((num_spots_batch_temp+1):sum(num_spots_batch),1) = merge_unique_variant_guides_full',num2str(j),'.variant_index_end_var;']);
                    eval(['batch_unique_variant_guides_full',num2str(j),'.chr_id_var((num_spots_batch_temp+1):sum(num_spots_batch),1) = merge_unique_variant_guides_full',num2str(j),'.chr_id_var;']);
                end
            elseif num_iterate == max(num_associated_variants_batch_temp)
                for h = 1:num_iterate;
                    eval(['batch_unique_variant_guides_full',num2str(h),'.bpstart_var = [batch_unique_variant_guides_full',num2str(h),'.bpstart_var;merge_unique_variant_guides_full',num2str(h),'.bpstart_var];']);
                    eval(['batch_unique_variant_guides_full',num2str(h),'.bpend_var = [batch_unique_variant_guides_full',num2str(h),'.bpend_var;merge_unique_variant_guides_full',num2str(h),'.bpend_var];']);
                    eval(['batch_unique_variant_guides_full',num2str(h),'.rs_id_var = [batch_unique_variant_guides_full',num2str(h),'.rs_id_var;merge_unique_variant_guides_full',num2str(h),'.rs_id_var];']);
                    eval(['batch_unique_variant_guides_full',num2str(h),'.ref_snp_var = [batch_unique_variant_guides_full',num2str(h),'.ref_snp_var;merge_unique_variant_guides_full',num2str(h),'.ref_snp_var];']);
                    eval(['batch_unique_variant_guides_full',num2str(h),'.alt_snp_var = [batch_unique_variant_guides_full',num2str(h),'.alt_snp_var;merge_unique_variant_guides_full',num2str(h),'.alt_snp_var];']);
                    eval(['batch_unique_variant_guides_full',num2str(h),'.identifier_var = [batch_unique_variant_guides_full',num2str(h),'.identifier_var;merge_unique_variant_guides_full',num2str(h),'.identifier_var];']);
                    eval(['batch_unique_variant_guides_full',num2str(h),'.variant_index_start_var = [batch_unique_variant_guides_full',num2str(h),'.variant_index_start_var;merge_unique_variant_guides_full',num2str(h),'.variant_index_start_var];']);
                    eval(['batch_unique_variant_guides_full',num2str(h),'.variant_index_end_var = [batch_unique_variant_guides_full',num2str(h),'.variant_index_end_var;merge_unique_variant_guides_full',num2str(h),'.variant_index_end_var];']);
                    eval(['batch_unique_variant_guides_full',num2str(h),'.chr_id_var = [batch_unique_variant_guides_full',num2str(h),'.chr_id_var;merge_unique_variant_guides_full',num2str(h),'.chr_id_var];']);
                end
            end
        end
        
        %% Haplotype Analysis Output
        
        if perform_haplotype_analysis == 1 && exist('unique_haplotype_guides_full','var')==1 && length(unique_haplotype_guides_full.guides)>0
            merge_output_filenames_list_haplo_temp(1:length(unique_haplotype_guides_full.guides)) = cellstr(temp_filename);
            merge_unique_haplotype_guides_full.guide_ID_number = [merge_unique_haplotype_guides_full.guide_ID_number;unique_haplotype_guides_full.guide_ID_number];
            merge_unique_haplotype_guides_full.guides_only = [merge_unique_haplotype_guides_full.guides_only;unique_haplotype_guides_full.guides_only];
            merge_unique_haplotype_guides_full.pam = [merge_unique_haplotype_guides_full.pam;unique_haplotype_guides_full.pam];
            merge_unique_haplotype_guides_full.pam_only = [merge_unique_haplotype_guides_full.pam_only;unique_haplotype_guides_full.pam_only];
            merge_unique_haplotype_guides_full.guides = [merge_unique_haplotype_guides_full.guides;unique_haplotype_guides_full.guides];
            merge_unique_haplotype_guides_full.chr_id = [merge_unique_haplotype_guides_full.chr_id;unique_haplotype_guides_full.chr_id];
            merge_unique_haplotype_guides_full.strand = [merge_unique_haplotype_guides_full.strand;unique_haplotype_guides_full.strand];
            merge_unique_haplotype_guides_full.count = [merge_unique_haplotype_guides_full.count;unique_haplotype_guides_full.count];
            merge_unique_haplotype_guides_full.frequency = [merge_unique_haplotype_guides_full.frequency;unique_haplotype_guides_full.frequency];
            merge_unique_haplotype_guides_full.filename = [merge_unique_haplotype_guides_full.filename;transpose(merge_output_filenames_list_haplo_temp)];
            
            merge_unique_haplotype_guides_full.pam_create = [merge_unique_haplotype_guides_full.pam_create;unique_haplotype_guides_full.pam_create];
            if exist('original_guides_full','var')==1
                merge_unique_haplotype_guides_full.n_mm_wt = [merge_unique_haplotype_guides_full.n_mm_wt;unique_haplotype_guides_full.n_mm_wt];
                merge_unique_haplotype_guides_full.guides_wt = [merge_unique_haplotype_guides_full.guides_wt;unique_haplotype_guides_full.guides_wt];
                merge_unique_haplotype_guides_full.guides_only_wt = [merge_unique_haplotype_guides_full.guides_only_wt;unique_haplotype_guides_full.guides_only_wt];
                merge_unique_haplotype_guides_full.guide_ID_number_wt = [merge_unique_haplotype_guides_full.guide_ID_number_wt;unique_haplotype_guides_full.guide_ID_number_wt];
                merge_unique_haplotype_guides_full.chr_id_wt = [merge_unique_haplotype_guides_full.chr_id_wt;unique_haplotype_guides_full.chr_id_wt];
                merge_unique_haplotype_guides_full.pam_only_wt = [merge_unique_haplotype_guides_full.pam_only_wt;unique_haplotype_guides_full.pam_only_wt];
                merge_unique_haplotype_guides_full.strand_wt = [merge_unique_haplotype_guides_full.strand_wt;unique_haplotype_guides_full.strand_wt];
                merge_unique_haplotype_guides_full.dsb_index_wt = [merge_unique_haplotype_guides_full.dsb_index_wt;unique_haplotype_guides_full.dsb_index_wt];
                merge_unique_haplotype_guides_full.dsb_coord_wt = [merge_unique_haplotype_guides_full.dsb_coord_wt;unique_haplotype_guides_full.dsb_coord_wt];
                merge_unique_haplotype_guides_full.guide_index_start_wt = [merge_unique_haplotype_guides_full.guide_index_start_wt;unique_haplotype_guides_full.guide_index_start_wt];
                merge_unique_haplotype_guides_full.guide_index_end_wt = [merge_unique_haplotype_guides_full.guide_index_end_wt;unique_haplotype_guides_full.guide_index_end_wt];
                merge_unique_haplotype_guides_full.guide_coord_start_wt = [merge_unique_haplotype_guides_full.guide_coord_start_wt;unique_haplotype_guides_full.guide_coord_start_wt];
                merge_unique_haplotype_guides_full.guide_coord_end_wt = [merge_unique_haplotype_guides_full.guide_coord_end_wt;unique_haplotype_guides_full.guide_coord_end_wt];
            end
            num_iterate = num_associated_variants(b,1);
            if max(num_associated_variants_temp) > num_iterate
                num_iterate = max(num_associated_variants_temp);
            end
            for h = 1:num_iterate;
                if eval(['exist([''merge_unique_haplotype_guides_full',num2str(h),'''],''var'')==0'])
                    
                    eval(['merge_unique_haplotype_guides_full',num2str(h),'.bpstart_var = [];']);
                    eval(['merge_unique_haplotype_guides_full',num2str(h),'.bpend_var = [];']);
                    eval(['merge_unique_haplotype_guides_full',num2str(h),'.rs_id_var = [];']);
                    eval(['merge_unique_haplotype_guides_full',num2str(h),'.ref_snp_var = [];']);
                    eval(['merge_unique_haplotype_guides_full',num2str(h),'.alt_snp_var = [];']);
                    eval(['merge_unique_haplotype_guides_full',num2str(h),'.identifier_var = [];']);
                    eval(['merge_unique_haplotype_guides_full',num2str(h),'.variant_index_start_var = [];']);
                    eval(['merge_unique_haplotype_guides_full',num2str(h),'.variant_index_end_var = [];']);
                    eval(['merge_unique_haplotype_guides_full',num2str(h),'.chr_id_var = [];']);
                end
            end
            num_iterate = num_associated_variants(b,1);
            if max(num_associated_variants_temp) > num_iterate
                num_iterate = max(num_associated_variants_temp);
                for j = (num_associated_variants(b,1)+1):max(num_associated_variants_temp)
                    if eval(['exist([''haplotype_associated_variants',num2str(j),'''],''var'')==0'])
                        eval(['haplotype_associated_variants',num2str(j),'.bpstart(1:num_spots(b,1),1) = {''''};']);
                        eval(['haplotype_associated_variants',num2str(j),'.bpend(1:num_spots(b,1),1) = {''''};']);
                        eval(['haplotype_associated_variants',num2str(j),'.rs_id(1:num_spots(b,1),1) = {''''};']);
                        eval(['haplotype_associated_variants',num2str(j),'.ref_snp(1:num_spots(b,1),1) = {''''};']);
                        eval(['haplotype_associated_variants',num2str(j),'.alt_snp(1:num_spots(b,1),1) = {''''};']);
                        eval(['haplotype_associated_variants',num2str(j),'.identifier(1:num_spots(b,1),1) = {''''};']);
                        eval(['haplotype_associated_variants',num2str(j),'.variant_index_start(1:num_spots(b,1),1) = {''''};']);
                        eval(['haplotype_associated_variants',num2str(j),'.variant_index_end(1:num_spots(b,1),1) = {''''};']);
                        eval(['haplotype_associated_variants',num2str(j),'.chr_id(1:num_spots(b,1),1) = {''''};']);
                    end
                end
            end
            if num_iterate > max(num_associated_variants_temp)
                for j = (max(num_associated_variants_temp)+1):num_associated_variants(b,1)
                    eval(['haplotype_associated_variants',num2str(j),'.bpstart((num_spots_sum(b-1)+1):num_spots_sum(b,1),1) = haplotype_associated_variants',num2str(h),'.bpstart;']);
                    eval(['haplotype_associated_variants',num2str(j),'.bpend((num_spots_sum(b-1)+1):num_spots_sum(b,1),1) = haplotype_associated_variants',num2str(h),'.bpend;']);
                    eval(['haplotype_associated_variants',num2str(j),'.rs_id((num_spots_sum(b-1)+1):num_spots_sum(b,1),1) = haplotype_associated_variants',num2str(h),'.rs_id;']);
                    eval(['haplotype_associated_variants',num2str(j),'.ref_snp((num_spots_sum(b-1)+1):num_spots_sum(b,1),1) = haplotype_associated_variants',num2str(h),'.ref_snp;']);
                    eval(['haplotype_associated_variants',num2str(j),'.alt_snp((num_spots_sum(b-1)+1):num_spots_sum(b,1),1) = haplotype_associated_variants',num2str(h),'.alt_snp;']);
                    eval(['haplotype_associated_variants',num2str(j),'.identifier((num_spots_sum(b-1)+1):num_spots_sum(b,1),1) = haplotype_associated_variants',num2str(h),'.identifier;']);
                    eval(['haplotype_associated_variants',num2str(j),'.variant_index_start((num_spots_sum(b-1)+1):num_spots_sum(b,1),1) = haplotype_associated_variants',num2str(h),'.variant_index_start;']);
                    eval(['haplotype_associated_variants',num2str(j),'.variant_index_end((num_spots_sum(b-1)+1):num_spots_sum(b,1),1) = haplotype_associated_variants',num2str(h),'.variant_index_end;']);
                    eval(['haplotype_associated_variants',num2str(j),'.chr_id((num_spots_sum(b-1)+1):num_spots_sum(b,1),1) = haplotype_associated_variants',num2str(h),'.chr_id;']);
                    
                    eval(['haplotype_associated_variants',num2str(j),'.bpstart(1:num_spots(b,1),1) = {''''};;']);
                    eval(['haplotype_associated_variants',num2str(j),'.bpend(1:num_spots(b,1),1) = {''''};;']);
                    eval(['haplotype_associated_variants',num2str(j),'.rs_id(1:num_spots(b,1),1) = {''''};;']);
                    eval(['haplotype_associated_variants',num2str(j),'.ref_snp(1:num_spots(b,1),1) = {''''};;']);
                    eval(['haplotype_associated_variants',num2str(j),'.alt_snp(1:num_spots(b,1),1) = {''''};;']);
                    eval(['haplotype_associated_variants',num2str(j),'.identifier(1:num_spots(b,1),1) = {''''};;']);
                    eval(['haplotype_associated_variants',num2str(j),'.variant_index_start(1:num_spots(b,1),1) = {''''};;']);
                    eval(['haplotype_associated_variants',num2str(j),'.variant_index_end(1:num_spots(b,1),1) = {''''};;']);
                    eval(['haplotype_associated_variants',num2str(j),'.chr_id(1:num_spots(b,1),1) = {''''};;']);
                end
            end
            
            for h = 1:num_iterate;
                eval(['merge_unique_haplotype_guides_full',num2str(h),'.bpstart_var = [merge_unique_haplotype_guides_full',num2str(h),'.bpstart_var;haplotype_associated_variants',num2str(h),'.bpstart];']);
                eval(['merge_unique_haplotype_guides_full',num2str(h),'.bpend_var = [merge_unique_haplotype_guides_full',num2str(h),'.bpend_var;haplotype_associated_variants',num2str(h),'.bpend];']);
                eval(['merge_unique_haplotype_guides_full',num2str(h),'.rs_id_var = [merge_unique_haplotype_guides_full',num2str(h),'.rs_id_var;haplotype_associated_variants',num2str(h),'.rs_id];']);
                eval(['merge_unique_haplotype_guides_full',num2str(h),'.ref_snp_var = [merge_unique_haplotype_guides_full',num2str(h),'.ref_snp_var;haplotype_associated_variants',num2str(h),'.ref_snp];']);
                eval(['merge_unique_haplotype_guides_full',num2str(h),'.alt_snp_var = [merge_unique_haplotype_guides_full',num2str(h),'.alt_snp_var;haplotype_associated_variants',num2str(h),'.alt_snp];']);
                eval(['merge_unique_haplotype_guides_full',num2str(h),'.identifier_var = [merge_unique_haplotype_guides_full',num2str(h),'.identifier_var;haplotype_associated_variants',num2str(h),'.identifier];']);
                eval(['merge_unique_haplotype_guides_full',num2str(h),'.variant_index_start_var = [merge_unique_haplotype_guides_full',num2str(h),'.variant_index_start_var;haplotype_associated_variants',num2str(h),'.variant_index_start];']);
                eval(['merge_unique_haplotype_guides_full',num2str(h),'.variant_index_end_var = [merge_unique_haplotype_guides_full',num2str(h),'.variant_index_end_var;haplotype_associated_variants',num2str(h),'.variant_index_end];']);
                eval(['merge_unique_haplotype_guides_full',num2str(h),'.chr_id_var = [merge_unique_haplotype_guides_full',num2str(h),'.chr_id_var;haplotype_associated_variants',num2str(h),'.chr_id];']);
            end
            
            
            if b == length(PAM_sequence_list)
                num_spots_batch(a,1) = sum(num_spots);
                num_spots_batch_temp = num_spots_batch;
                num_spots_batch_temp(a,:) = [];
                if length(num_spots_batch) == 1
                    num_spots_batch_temp = 0;
                elseif length(num_spots_batch)>1
                    num_spots_batch_temp = sum(num_spots_batch_temp);
                end
            end
        end
        if end_number > starting_number && perform_batch_output == 1 && perform_haplotype_analysis == 1 && exist('merge_unique_haplotype_guides_full','var')==1 && b == length(PAM_sequence_list)
            batch_output_filenames_list_haplo_temp(1:length(merge_unique_haplotype_guides_full.guides)) = cellstr(temp_filename);
            batch_unique_haplotype_guides_full.guide_ID_number = [batch_unique_haplotype_guides_full.guide_ID_number;merge_unique_haplotype_guides_full.guide_ID_number];
            batch_unique_haplotype_guides_full.guides_only = [batch_unique_haplotype_guides_full.guides_only;merge_unique_haplotype_guides_full.guides_only];
            batch_unique_haplotype_guides_full.pam = [batch_unique_haplotype_guides_full.pam;merge_unique_haplotype_guides_full.pam];
            batch_unique_haplotype_guides_full.pam_only = [batch_unique_haplotype_guides_full.pam_only;merge_unique_haplotype_guides_full.pam_only];
            batch_unique_haplotype_guides_full.guides = [batch_unique_haplotype_guides_full.guides;merge_unique_haplotype_guides_full.guides];
            batch_unique_haplotype_guides_full.chr_id = [batch_unique_haplotype_guides_full.chr_id;merge_unique_haplotype_guides_full.chr_id];
            batch_unique_haplotype_guides_full.strand = [batch_unique_haplotype_guides_full.strand;merge_unique_haplotype_guides_full.strand];
            batch_unique_haplotype_guides_full.count = [batch_unique_haplotype_guides_full.count;merge_unique_haplotype_guides_full.count];
            batch_unique_haplotype_guides_full.frequency = [batch_unique_haplotype_guides_full.frequency;merge_unique_haplotype_guides_full.frequency];
            batch_unique_haplotype_guides_full.filename = [batch_unique_haplotype_guides_full.filename;transpose(batch_output_filenames_list_haplo_temp)];
            
            batch_unique_haplotype_guides_full.pam_create = [batch_unique_haplotype_guides_full.pam_create;merge_unique_haplotype_guides_full.pam_create];
            batch_unique_haplotype_guides_full.n_mm_wt = [batch_unique_haplotype_guides_full.n_mm_wt;merge_unique_haplotype_guides_full.n_mm_wt];
            batch_unique_haplotype_guides_full.guides_wt = [batch_unique_haplotype_guides_full.guides_wt;merge_unique_haplotype_guides_full.guides_wt];
            batch_unique_haplotype_guides_full.guides_only_wt = [batch_unique_haplotype_guides_full.guides_only_wt;merge_unique_haplotype_guides_full.guides_only_wt];
            batch_unique_haplotype_guides_full.guide_ID_number_wt = [batch_unique_haplotype_guides_full.guide_ID_number_wt;merge_unique_haplotype_guides_full.guide_ID_number_wt];
            batch_unique_haplotype_guides_full.chr_id_wt = [batch_unique_haplotype_guides_full.chr_id_wt;merge_unique_haplotype_guides_full.chr_id_wt];
            batch_unique_haplotype_guides_full.pam_only_wt = [batch_unique_haplotype_guides_full.pam_only_wt;merge_unique_haplotype_guides_full.pam_only_wt];
            batch_unique_haplotype_guides_full.strand_wt = [batch_unique_haplotype_guides_full.strand_wt;merge_unique_haplotype_guides_full.strand_wt];
            batch_unique_haplotype_guides_full.dsb_index_wt = [batch_unique_haplotype_guides_full.dsb_index_wt;merge_unique_haplotype_guides_full.dsb_index_wt];
            batch_unique_haplotype_guides_full.dsb_coord_wt = [batch_unique_haplotype_guides_full.dsb_coord_wt;merge_unique_haplotype_guides_full.dsb_coord_wt];
            batch_unique_haplotype_guides_full.guide_index_start_wt = [batch_unique_haplotype_guides_full.guide_index_start_wt;merge_unique_haplotype_guides_full.guide_index_start_wt];
            batch_unique_haplotype_guides_full.guide_index_end_wt = [batch_unique_haplotype_guides_full.guide_index_end_wt;merge_unique_haplotype_guides_full.guide_index_end_wt];
            batch_unique_haplotype_guides_full.guide_coord_start_wt = [batch_unique_haplotype_guides_full.guide_coord_start_wt;merge_unique_haplotype_guides_full.guide_coord_start_wt];
            batch_unique_haplotype_guides_full.guide_coord_end_wt = [batch_unique_haplotype_guides_full.guide_coord_end_wt;merge_unique_haplotype_guides_full.guide_coord_end_wt];
            
            num_iterate = num_associated_variants_batch(a,1);
            num_associated_variants_batch_temp = num_associated_variants_batch;
            if a > 1
                num_associated_variants_batch_temp(a,:) = [];
            end
            
            for h = 1:num_iterate;
                if eval(['exist([''batch_unique_haplotype_guides_full',num2str(h),'''],''var'')==0'])
                    eval(['batch_unique_haplotype_guides_full',num2str(h),'.bpstart_var = {};']);
                    eval(['batch_unique_haplotype_guides_full',num2str(h),'.bpend_var = {};']);
                    eval(['batch_unique_haplotype_guides_full',num2str(h),'.rs_id_var = {};']);
                    eval(['batch_unique_haplotype_guides_full',num2str(h),'.ref_snp_var = {};']);
                    eval(['batch_unique_haplotype_guides_full',num2str(h),'.alt_snp_var = {};']);
                    eval(['batch_unique_haplotype_guides_full',num2str(h),'.identifier_var = {};']);
                    eval(['batch_unique_haplotype_guides_full',num2str(h),'.variant_index_start_var = {};']);
                    eval(['batch_unique_haplotype_guides_full',num2str(h),'.variant_index_end_var = {};']);
                    eval(['batch_unique_haplotype_guides_full',num2str(h),'.chr_id_var = {};']);
                end
            end
            
            num_iterate = num_associated_variants_batch(a,1);
            if max(num_associated_variants_batch_temp) > num_iterate
                for j = (num_associated_variants_batch(a,1)+1):max(num_associated_variants_batch_temp)
                    eval(['merge_unique_haplotype_guides_full',num2str(j),'.bpstart_var(1:num_spots_batch(a,1),1) = {''''};']);
                    eval(['merge_unique_haplotype_guides_full',num2str(j),'.bpend_var(1:num_spots_batch(a,1),1) = {''''};']);
                    eval(['merge_unique_haplotype_guides_full',num2str(j),'.rs_id_var(1:num_spots_batch(a,1),1) = {''''};']);
                    eval(['merge_unique_haplotype_guides_full',num2str(j),'.ref_snp_var(1:num_spots_batch(a,1),1) = {''''};']);
                    eval(['merge_unique_haplotype_guides_full',num2str(j),'.alt_snp_var(1:num_spots_batch(a,1),1) = {''''};']);
                    eval(['merge_unique_haplotype_guides_full',num2str(j),'.identifier_var(1:num_spots_batch(a,1),1) = {''''};']);
                    eval(['merge_unique_haplotype_guides_full',num2str(j),'.variant_index_start_var(1:num_spots_batch(a,1),1) = {''''};']);
                    eval(['merge_unique_haplotype_guides_full',num2str(j),'.variant_index_end_var(1:num_spots_batch(a,1),1) = {''''};']);
                    eval(['merge_unique_haplotype_guides_full',num2str(j),'.chr_id_var(1:num_spots_batch(a,1),1) = {''''};']);
                end
                for h = 1:max(num_associated_variants_batch_temp);
                    eval(['batch_unique_haplotype_guides_full',num2str(h),'.bpstart_var = [batch_unique_haplotype_guides_full',num2str(h),'.bpstart_var;merge_unique_haplotype_guides_full',num2str(h),'.bpstart_var];']);
                    eval(['batch_unique_haplotype_guides_full',num2str(h),'.bpend_var = [batch_unique_haplotype_guides_full',num2str(h),'.bpend_var;merge_unique_haplotype_guides_full',num2str(h),'.bpend_var];']);
                    eval(['batch_unique_haplotype_guides_full',num2str(h),'.rs_id_var = [batch_unique_haplotype_guides_full',num2str(h),'.rs_id_var;merge_unique_haplotype_guides_full',num2str(h),'.rs_id_var];']);
                    eval(['batch_unique_haplotype_guides_full',num2str(h),'.ref_snp_var = [batch_unique_haplotype_guides_full',num2str(h),'.ref_snp_var;merge_unique_haplotype_guides_full',num2str(h),'.ref_snp_var];']);
                    eval(['batch_unique_haplotype_guides_full',num2str(h),'.alt_snp_var = [batch_unique_haplotype_guides_full',num2str(h),'.alt_snp_var;merge_unique_haplotype_guides_full',num2str(h),'.alt_snp_var];']);
                    eval(['batch_unique_haplotype_guides_full',num2str(h),'.identifier_var = [batch_unique_haplotype_guides_full',num2str(h),'.identifier_var;merge_unique_haplotype_guides_full',num2str(h),'.identifier_var];']);
                    eval(['batch_unique_haplotype_guides_full',num2str(h),'.variant_index_start_var = [batch_unique_haplotype_guides_full',num2str(h),'.variant_index_start_var;merge_unique_haplotype_guides_full',num2str(h),'.variant_index_start_var];']);
                    eval(['batch_unique_haplotype_guides_full',num2str(h),'.variant_index_end_var = [batch_unique_haplotype_guides_full',num2str(h),'.variant_index_end_var;merge_unique_haplotype_guides_full',num2str(h),'.variant_index_end_var];']);
                    eval(['batch_unique_haplotype_guides_full',num2str(h),'.chr_id_var = [batch_unique_haplotype_guides_full',num2str(h),'.chr_id_var;merge_unique_haplotype_guides_full',num2str(h),'.chr_id_var];']);
                end
            elseif num_iterate > max(num_associated_variants_batch_temp)
                for j = 1:max(num_associated_variants_batch(a,1))
                    eval(['batch_unique_haplotype_guides_full',num2str(j),'.bpstart_var((num_spots_batch_temp+1):sum(num_spots_batch),1) = merge_unique_haplotype_guides_full',num2str(j),'.bpstart_var;']);
                    eval(['batch_unique_haplotype_guides_full',num2str(j),'.bpend_var((num_spots_batch_temp+1):sum(num_spots_batch),1) = merge_unique_haplotype_guides_full',num2str(j),'.bpend_var;']);
                    eval(['batch_unique_haplotype_guides_full',num2str(j),'.rs_id_var((num_spots_batch_temp+1):sum(num_spots_batch),1) = merge_unique_haplotype_guides_full',num2str(j),'.rs_id_var;']);
                    eval(['batch_unique_haplotype_guides_full',num2str(j),'.ref_snp_var((num_spots_batch_temp+1):sum(num_spots_batch),1) = merge_unique_haplotype_guides_full',num2str(j),'.ref_snp_var;']);
                    eval(['batch_unique_haplotype_guides_full',num2str(j),'.alt_snp_var((num_spots_batch_temp+1):sum(num_spots_batch),1) = merge_unique_haplotype_guides_full',num2str(j),'.alt_snp_var;']);
                    eval(['batch_unique_haplotype_guides_full',num2str(j),'.identifier_var((num_spots_batch_temp+1):sum(num_spots_batch),1) = merge_unique_haplotype_guides_full',num2str(j),'.identifier_var;']);
                    eval(['batch_unique_haplotype_guides_full',num2str(j),'.variant_index_start_var((num_spots_batch_temp+1):sum(num_spots_batch),1) = merge_unique_haplotype_guides_full',num2str(j),'.variant_index_start_var;']);
                    eval(['batch_unique_haplotype_guides_full',num2str(j),'.variant_index_end_var((num_spots_batch_temp+1):sum(num_spots_batch),1) = merge_unique_haplotype_guides_full',num2str(j),'.variant_index_end_var;']);
                    eval(['batch_unique_haplotype_guides_full',num2str(j),'.chr_id_var((num_spots_batch_temp+1):sum(num_spots_batch),1) = merge_unique_haplotype_guides_full',num2str(j),'.chr_id_var;']);
                end
            elseif num_iterate == max(num_associated_variants_batch_temp)
                for h = 1:num_iterate;
                    eval(['batch_unique_haplotype_guides_full',num2str(h),'.bpstart_var = [batch_unique_haplotype_guides_full',num2str(h),'.bpstart_var;merge_unique_haplotype_guides_full',num2str(h),'.bpstart_var];']);
                    eval(['batch_unique_haplotype_guides_full',num2str(h),'.bpend_var = [batch_unique_haplotype_guides_full',num2str(h),'.bpend_var;merge_unique_haplotype_guides_full',num2str(h),'.bpend_var];']);
                    eval(['batch_unique_haplotype_guides_full',num2str(h),'.rs_id_var = [batch_unique_haplotype_guides_full',num2str(h),'.rs_id_var;merge_unique_haplotype_guides_full',num2str(h),'.rs_id_var];']);
                    eval(['batch_unique_haplotype_guides_full',num2str(h),'.ref_snp_var = [batch_unique_haplotype_guides_full',num2str(h),'.ref_snp_var;merge_unique_haplotype_guides_full',num2str(h),'.ref_snp_var];']);
                    eval(['batch_unique_haplotype_guides_full',num2str(h),'.alt_snp_var = [batch_unique_haplotype_guides_full',num2str(h),'.alt_snp_var;merge_unique_haplotype_guides_full',num2str(h),'.alt_snp_var];']);
                    eval(['batch_unique_haplotype_guides_full',num2str(h),'.identifier_var = [batch_unique_haplotype_guides_full',num2str(h),'.identifier_var;merge_unique_haplotype_guides_full',num2str(h),'.identifier_var];']);
                    eval(['batch_unique_haplotype_guides_full',num2str(h),'.variant_index_start_var = [batch_unique_haplotype_guides_full',num2str(h),'.variant_index_start_var;merge_unique_haplotype_guides_full',num2str(h),'.variant_index_start_var];']);
                    eval(['batch_unique_haplotype_guides_full',num2str(h),'.variant_index_end_var = [batch_unique_haplotype_guides_full',num2str(h),'.variant_index_end_var;merge_unique_haplotype_guides_full',num2str(h),'.variant_index_end_var];']);
                    eval(['batch_unique_haplotype_guides_full',num2str(h),'.chr_id_var = [batch_unique_haplotype_guides_full',num2str(h),'.chr_id_var;merge_unique_haplotype_guides_full',num2str(h),'.chr_id_var];']);
                end
            end
        end
        alphabet_letters1 = cellstr(['A';'B';'C';'D';'E';'F';'G';'H';'I';'J';'K';'L';'M';'N';'O';'P';'Q';'R';'S';'T';'U';'V';'W';'X';'Y';'Z']);
        alphabet_letters2 = cellstr(['AA';'AB';'AC';'AD';'AE';'AF';'AG';'AH';'AI';'AJ';'AK';'AL';'AM';'AN';'AO';'AP';'AQ';'AR';'AS';'AT';'AU';'AV';'AW';'AX';'AY';'AZ']);
        alphabet_letters3 = cellstr(['BA';'BB';'BC';'BD';'BE';'BF';'BG';'BH';'BI';'BJ';'BK';'BL';'BM';'BN';'BO';'BP';'BQ';'BR';'BS';'BT';'BU';'BV';'BW';'BX';'BY';'BZ']);
        alphabet_letters4 = cellstr(['CA';'CB';'CC';'CD';'CE';'CF';'CG';'CH';'CI';'CJ';'CK';'CL';'CM';'CN';'CO';'CP';'CQ';'CR';'CS';'CT';'CU';'CV';'CW';'CX';'CY';'CZ']);
        alphabet_letters = vertcat(alphabet_letters1,alphabet_letters2,alphabet_letters3,alphabet_letters4);
        
        for v = 1:length(alphabet_letters);
            alphabet1(v,1) = cellstr(strjoin([alphabet_letters(v),'1'],''));
            alphabet2(v,1) = cellstr(strjoin([alphabet_letters(v),'2'],''));
        end
        
        if end_number == starting_number && output_to_excel == 1 && b==length(PAM_sequence_list) || individual_output == 1 && output_to_excel == 1 && b==length(PAM_sequence_list)
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'ID'}, 'No Variant Analysis', char(alphabet1(1)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_original_guides_full.guide_ID_number, 'No Variant Analysis', char(alphabet2(1)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'PAM'}, 'No Variant Analysis', char(alphabet1(2)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_original_guides_full.pam, 'No Variant Analysis', char(alphabet2(2)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'sgRNA Sequence'}, 'No Variant Analysis', char(alphabet1(3)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_original_guides_full.guides_only, 'No Variant Analysis', char(alphabet2(3)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'PAM Only'}, 'No Variant Analysis', char(alphabet1(4)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_original_guides_full.pam_only, 'No Variant Analysis', char(alphabet2(4)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'sgRNA+PAM'}, 'No Variant Analysis', char(alphabet1(5)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_original_guides_full.guides, 'No Variant Analysis', char(alphabet2(5)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'Strand'}, 'No Variant Analysis', char(alphabet1(6)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_original_guides_full.strand, 'No Variant Analysis', char(alphabet2(6)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'Chromosome'}, 'No Variant Analysis', char(alphabet1(7)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_original_guides_full.chr_id, 'No Variant Analysis', char(alphabet2(7)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], temp_coord_string_start, 'No Variant Analysis', char(alphabet1(8)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_original_guides_full.guide_coord_start, 'No Variant Analysis', char(alphabet2(8)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], temp_coord_string_end, 'No Variant Analysis', char(alphabet1(9)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_original_guides_full.guide_coord_end, 'No Variant Analysis', char(alphabet2(9)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], temp_coord_string_dsb, 'No Variant Analysis', char(alphabet1(10)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_original_guides_full.dsb_coord, 'No Variant Analysis', char(alphabet2(10)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'Guide Start (Index)'}, 'No Variant Analysis', char(alphabet1(11)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_original_guides_full.guide_index_start, 'No Variant Analysis', char(alphabet2(11)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'Guide End (Index)'}, 'No Variant Analysis', char(alphabet1(12)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_original_guides_full.guide_index_end, 'No Variant Analysis', char(alphabet2(12)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'Double Strand Break (Index)'}, 'No Variant Analysis', char(alphabet1(13)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_original_guides_full.dsb_index, 'No Variant Analysis', char(alphabet2(13)));
            if multiple_match_analysis == 1
                xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'Predicted # of DSBs in Sequence'}, 'No Variant Analysis', char(alphabet1(14)));
                xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_original_guides_full.mult_match_count, 'No Variant Analysis', char(alphabet2(14)));
            end
        end
        
        if multiple_match_analysis == 1 && end_number == starting_number && output_to_excel == 1 && exist('multiple_match','var')==1 && b==length(PAM_sequence_list) || multiple_match_analysis == 1 && individual_output == 1 && output_to_excel == 1 && exist('multiple_match','var')>0 && b==length(PAM_sequence_list)
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'ID'}, 'Multiple Match Analysis', char(alphabet1(1)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_multiple_match.guide_ID_number, 'Multiple Match Analysis', char(alphabet2(1)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'PAM'}, 'Multiple Match Analysis', char(alphabet1(2)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_multiple_match.pam, 'Multiple Match Analysis', char(alphabet2(2)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'sgRNA Sequence'}, 'Multiple Match Analysis', char(alphabet1(3)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_multiple_match.guides_only, 'Multiple Match Analysis', char(alphabet2(3)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'PAM Only'}, 'Multiple Match Analysis', char(alphabet1(4)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_multiple_match.pam_only, 'Multiple Match Analysis', char(alphabet2(4)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'sgRNA+PAM'}, 'Multiple Match Analysis', char(alphabet1(5)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_multiple_match.guides, 'Multiple Match Analysis', char(alphabet2(5)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'Strand'}, 'Multiple Match Analysis', char(alphabet1(6)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_multiple_match.strand, 'Multiple Match Analysis', char(alphabet2(6)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'Chromosome'}, 'Multiple Match Analysis', char(alphabet1(7)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_multiple_match.chr_id, 'Multiple Match Analysis', char(alphabet2(7)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], temp_coord_string_start, 'Multiple Match Analysis', char(alphabet1(8)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_multiple_match.guide_coord_start, 'Multiple Match Analysis', char(alphabet2(8)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], temp_coord_string_end, 'Multiple Match Analysis', char(alphabet1(9)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_multiple_match.guide_coord_end, 'Multiple Match Analysis', char(alphabet2(9)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], temp_coord_string_dsb, 'Multiple Match Analysis', char(alphabet1(10)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_multiple_match.dsb_coord, 'Multiple Match Analysis', char(alphabet2(10)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'Guide Start (Index)'}, 'Multiple Match Analysis', char(alphabet1(11)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_multiple_match.guide_index_start, 'Multiple Match Analysis', char(alphabet2(11)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'Guide End (Index)'}, 'Multiple Match Analysis', char(alphabet1(12)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_multiple_match.guide_index_end, 'Multiple Match Analysis', char(alphabet2(12)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'Double Strand Break (Index)'}, 'Multiple Match Analysis', char(alphabet1(13)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_multiple_match.dsb_index, 'Multiple Match Analysis', char(alphabet2(13)));
        end
        
        if perform_variant_analysis == 1 && end_number == starting_number && output_to_excel == 1 && b==length(PAM_sequence_list) || perform_variant_analysis == 1 && individual_output == 1 && output_to_excel == 1 && b==length(PAM_sequence_list) || perform_wgs_analysis == 1 && end_number == starting_number && output_to_excel == 1 && b==length(PAM_sequence_list) || perform_wgs_analysis == 1 && individual_output == 1 && output_to_excel == 1 && b==length(PAM_sequence_list)
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'ID'}, 'Variant Analysis sgRNA', char(alphabet1(1)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_unique_variant_guides_full.guide_ID_number, 'Variant Analysis sgRNA', char(alphabet2(1)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'PAM'}, 'Variant Analysis sgRNA', char(alphabet1(2)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_unique_variant_guides_full.pam, 'Variant Analysis sgRNA', char(alphabet2(2)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'sgRNA Sequence'}, 'Variant Analysis sgRNA', char(alphabet1(3)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_unique_variant_guides_full.guides_only, 'Variant Analysis sgRNA', char(alphabet2(3)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'PAM'}, 'Variant Analysis sgRNA', char(alphabet1(4)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_unique_variant_guides_full.pam_only, 'Variant Analysis sgRNA', char(alphabet2(4)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'sgRNA+PAM'}, 'Variant Analysis sgRNA', char(alphabet1(5)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_unique_variant_guides_full.guides, 'Variant Analysis sgRNA', char(alphabet2(5)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'Strand'}, 'Variant Analysis sgRNA', char(alphabet1(6)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_unique_variant_guides_full.strand, 'Variant Analysis sgRNA', char(alphabet2(6)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'Chromosome'}, 'Variant Analysis sgRNA', char(alphabet1(7)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_unique_variant_guides_full.chr_id, 'Variant Analysis sgRNA', char(alphabet2(7)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'sgRNA Classification'}, 'Variant Analysis sgRNA', char(alphabet1(8)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_unique_variant_guides_full.pam_create, 'Variant Analysis sgRNA', char(alphabet2(8)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'N-W Alignment Score'}, 'Variant Analysis sgRNA', char(alphabet1(9)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_unique_variant_guides_full.n_mm_wt, 'Variant Analysis sgRNA', char(alphabet2(9)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'ID'}, 'Variant Analysis sgRNA', char(alphabet1(10)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_unique_variant_guides_full.guide_ID_number_wt, 'Variant Analysis sgRNA', char(alphabet2(10)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'Sequence of Associated No Variant sgRNA'}, 'Variant Analysis sgRNA', char(alphabet1(11)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_unique_variant_guides_full.guides_only_wt, 'Variant Analysis sgRNA', char(alphabet2(11)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'PAM'}, 'Variant Analysis sgRNA', char(alphabet1(12)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_unique_variant_guides_full.pam_only_wt, 'Variant Analysis sgRNA', char(alphabet2(12)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'sgRNA+PAM'}, 'Variant Analysis sgRNA', char(alphabet1(13)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_unique_variant_guides_full.guides_wt, 'Variant Analysis sgRNA', char(alphabet2(13)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'Strand'}, 'Variant Analysis sgRNA', char(alphabet1(14)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_unique_variant_guides_full.strand_wt, 'Variant Analysis sgRNA', char(alphabet2(14)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'Chromosome'}, 'Variant Analysis sgRNA', char(alphabet1(15)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_unique_variant_guides_full.chr_id_wt, 'Variant Analysis sgRNA', char(alphabet2(15)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], temp_coord_string_start, 'Variant Analysis sgRNA', char(alphabet1(16)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_unique_variant_guides_full.guide_coord_start_wt, 'Variant Analysis sgRNA', char(alphabet2(16)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], temp_coord_string_end, 'Variant Analysis sgRNA', char(alphabet1(17)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_unique_variant_guides_full.guide_coord_end_wt, 'Variant Analysis sgRNA', char(alphabet2(17)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], temp_coord_string_dsb, 'Variant Analysis sgRNA', char(alphabet1(18)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_unique_variant_guides_full.dsb_coord_wt, 'Variant Analysis sgRNA', char(alphabet2(18)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'Guide Start (Index)'}, 'Variant Analysis sgRNA', char(alphabet1(19)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_unique_variant_guides_full.guide_index_start_wt, 'Variant Analysis sgRNA', char(alphabet2(19)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'Guide End (Index)'}, 'Variant Analysis sgRNA', char(alphabet1(20)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_unique_variant_guides_full.guide_index_end_wt, 'Variant Analysis sgRNA', char(alphabet2(20)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'Double Strand Break (Index)'}, 'Variant Analysis sgRNA', char(alphabet1(21)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_unique_variant_guides_full.dsb_index_wt, 'Variant Analysis sgRNA', char(alphabet2(21)));
            
            num_iterate = num_associated_variants_batch(a,1);
            
            for h = 1:num_iterate;
                xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'RS ID'}, 'Variant Analysis sgRNA', char(alphabet1(24+9*(h-1))));
                xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], eval(['merge_unique_variant_guides_full',num2str(h),'.rs_id_var']), 'Variant Analysis sgRNA', char(alphabet2(24+9*(h-1))));
                xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'Ref SNP'}, 'Variant Analysis sgRNA', char(alphabet1(25+9*(h-1))));
                xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], eval(['merge_unique_variant_guides_full',num2str(h),'.ref_snp_var']), 'Variant Analysis sgRNA', char(alphabet2(25+9*(h-1))));
                xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'Alt SNP'}, 'Variant Analysis sgRNA', char(alphabet1(26+9*(h-1))));
                xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], eval(['merge_unique_variant_guides_full',num2str(h),'.alt_snp_var']), 'Variant Analysis sgRNA', char(alphabet2(26+9*(h-1))));
                xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'Variant Class'}, 'Variant Analysis sgRNA', char(alphabet1(27+9*(h-1))));
                xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], eval(['merge_unique_variant_guides_full',num2str(h),'.identifier_var']), 'Variant Analysis sgRNA', char(alphabet2(27+9*(h-1))));
                xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'Variant Chr'}, 'Variant Analysis sgRNA', char(alphabet1(28+9*(h-1))));
                xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], eval(['merge_unique_variant_guides_full',num2str(h),'.chr_id_var']), 'Variant Analysis sgRNA', char(alphabet2(28+9*(h-1))));
                xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], temp_coord_string_start_variant, 'Variant Analysis sgRNA', char(alphabet1(29+9*(h-1))));
                xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], eval(['merge_unique_variant_guides_full',num2str(h),'.bpstart_var']), 'Variant Analysis sgRNA', char(alphabet2(29+9*(h-1))));
                xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], temp_coord_string_end_variant, 'Variant Analysis sgRNA', char(alphabet1(30+9*(h-1))));
                xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], eval(['merge_unique_variant_guides_full',num2str(h),'.bpend_var']), 'Variant Analysis sgRNA', char(alphabet2(30+9*(h-1))));
                xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'Variant Start (Index)'}, 'Variant Analysis sgRNA', char(alphabet1(31+9*(h-1))));
                xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], eval(['merge_unique_variant_guides_full',num2str(h),'.variant_index_start_var']), 'Variant Analysis sgRNA', char(alphabet2(31+9*(h-1))));
                xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'Variant End (Index)'}, 'Variant Analysis sgRNA', char(alphabet1(32+9*(h-1))));
                xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], eval(['merge_unique_variant_guides_full',num2str(h),'.variant_index_end_var']), 'Variant Analysis sgRNA', char(alphabet2(32+9*(h-1))));
            end
        end
        
        if end_number == starting_number && output_to_excel == 0 && b==length(PAM_sequence_list) || individual_output == 1 && output_to_excel == 0 && b==length(PAM_sequence_list)
            array_cell_no_variant = {'ID' 'PAM' 'sgRNA_Sequence' 'PAM_Only' 'sgRNA_PAM' 'Strand' 'Chromosome' 'Guide_Start_Coordinate' 'Guide_End_Coordinate' 'DSB_Coordinate' 'Guide_Start_Index' 'Guide_End_Index' 'Double_Strand_Break_Index'};
            no_variant_output = table(merge_original_guides_full.guide_ID_number,merge_original_guides_full.pam,merge_original_guides_full.guides_only,merge_original_guides_full.pam_only,merge_original_guides_full.guides,merge_original_guides_full.strand,merge_original_guides_full.chr_id,merge_original_guides_full.guide_coord_start,merge_original_guides_full.guide_coord_end,merge_original_guides_full.dsb_coord,merge_original_guides_full.guide_index_start,merge_original_guides_full.guide_index_end,merge_original_guides_full.dsb_index,'VariableNames',array_cell_no_variant);
            writetable(no_variant_output,['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'_no_variant.txt'],'Delimiter','\t')
        end
        
        if multiple_match_analysis == 1 && end_number == starting_number && output_to_excel == 0 && b==length(PAM_sequence_list) && exist('multiple_match','var')==1 || multiple_match_analysis == 1 && individual_output == 1 && output_to_excel == 0 && exist('multiple_match','var')>0 && b==length(PAM_sequence_list)
            array_cell_multiple_match = {'ID' 'PAM' 'sgRNA_Sequence' 'PAM_Only' 'sgRNA_PAM' 'Strand' 'Chromosome' 'Guide_Start_Coordinate' 'Guide_End_Coordinate' 'DSB_Coordinate' 'Guide_Start_Index' 'Guide_End_Index' 'Double_Strand_Break_Index' 'Filename'};
            multiple_match_output = table(merge_multiple_match.guide_ID_number,merge_multiple_match.pam,merge_multiple_match.guides_only,merge_multiple_match.pam_only,merge_multiple_match.guides,merge_multiple_match.strand,merge_multiple_match.chr_id,merge_multiple_match.guide_coord_start,merge_multiple_match.guide_coord_end,merge_multiple_match.dsb_coord,merge_multiple_match.guide_index_start,merge_multiple_match.guide_index_end,merge_multiple_match.dsb_index,merge_multiple_match.filename,'VariableNames',array_cell_multiple_match);
            writetable(multiple_match_output,['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'_multiple_match.txt'],'Delimiter','\t')
        end
        
        if perform_haplotype_analysis == 1 && end_number == starting_number && b==length(PAM_sequence_list) && output_to_excel == 0 && exist('merge_unique_haplotype_guides_full','var')==1 || perform_haplotype_analysis == 1 && individual_output == 1 && output_to_excel == 0 && b==length(PAM_sequence_list) && exist('merge_unique_haplotype_guides_full','var')==1;
            array_cell_haplotype = {'ID' 'PAM' 'sgRNA_Sequence' 'sgRNA_PAM' 'Full_Sequence' 'Strand' 'Chromosome' 'Count' 'Frequency' 'Classification' 'N_W_Alignment_Score' 'No_Variant_ID' 'No_Variant_sgRNA_Sequence' 'No_Variant_PAM_Only' 'No_Variant_Full_Sequence' 'No_Variant_Strand' 'No_Variant_Chromosome' 'Guide_Start_Coordinate' 'Guide_End_Coordinate' 'DSB_Coordinate' 'Guide_Start_Index' 'Guide_End_Index' 'DSB_Index' 'Filename'};
            merge_haplotype_output = table(merge_unique_haplotype_guides_full.guide_ID_number,merge_unique_haplotype_guides_full.pam,merge_unique_haplotype_guides_full.guides_only,merge_unique_haplotype_guides_full.pam_only,merge_unique_haplotype_guides_full.guides,merge_unique_haplotype_guides_full.strand,merge_unique_haplotype_guides_full.chr_id,merge_unique_haplotype_guides_full.count,merge_unique_haplotype_guides_full.frequency,merge_unique_haplotype_guides_full.pam_create,merge_unique_haplotype_guides_full.n_mm_wt,merge_unique_haplotype_guides_full.guide_ID_number_wt,merge_unique_haplotype_guides_full.guides_only_wt,merge_unique_haplotype_guides_full.pam_only_wt,merge_unique_haplotype_guides_full.guides_wt,merge_unique_haplotype_guides_full.strand_wt,merge_unique_haplotype_guides_full.chr_id_wt,merge_unique_haplotype_guides_full.guide_coord_start_wt,merge_unique_haplotype_guides_full.guide_coord_end_wt,merge_unique_haplotype_guides_full.dsb_coord_wt,merge_unique_haplotype_guides_full.guide_index_start_wt,merge_unique_haplotype_guides_full.guide_index_end_wt,merge_unique_haplotype_guides_full.dsb_index_wt,merge_unique_haplotype_guides_full.filename,'VariableNames',array_cell_haplotype);
            
            num_iterate = num_associated_variants_batch(a,1);
            counter = 1;
            for h = 1:num_iterate
                array_cell_haplotype_variants(1,(counter:(counter+8))) = {['RS_ID_',num2str(h),''] ['Ref_SNP_',num2str(h),''] ['Alt_SNP_',num2str(h),''] ['Variant_Class_',num2str(h),''] ['Variant_Chromosome_',num2str(h),''] ['Variant_Start_Coordinate_',num2str(h),''] ['Variant_End_Coordinate_',num2str(h),''] ['Variant_Start_Index_',num2str(h),''] ['Variant_End_Index_',num2str(h),'']};
                counter = counter + 9;
            end
            counter = 1;
            for h = 1:num_iterate;
                eval(['haplotype_variants_output_temp',num2str(h),' = table(merge_unique_haplotype_guides_full',num2str(h),'.rs_id_var,merge_unique_haplotype_guides_full',num2str(h),'.ref_snp_var,merge_unique_haplotype_guides_full',num2str(h),'.alt_snp_var,merge_unique_haplotype_guides_full',num2str(h),'.identifier_var,merge_unique_haplotype_guides_full',num2str(h),'.chr_id_var,merge_unique_haplotype_guides_full',num2str(h),'.bpstart_var,merge_unique_haplotype_guides_full',num2str(h),'.bpend_var,merge_unique_haplotype_guides_full',num2str(h),'.variant_index_start_var,merge_unique_haplotype_guides_full',num2str(h),'.variant_index_end_var,''VariableNames'',array_cell_haplotype_variants(1,(counter:counter+8)));']);
                counter = counter + 9;
            end
            
            for h = 1:num_iterate;
                if h == 1;
                    eval(['haplotype_analysis_output',num2str(h),' = horzcat(merge_haplotype_output,haplotype_variants_output_temp',num2str(h),');']);
                elseif h == 2
                    eval(['haplotype_analysis_output',num2str(h),' = horzcat(haplotype_analysis_output',num2str(h-1),',haplotype_variants_output_temp',num2str(h),');']);
                elseif h>2
                    eval(['haplotype_analysis_output',num2str(h),' = horzcat(haplotype_analysis_output',num2str(h-1),',haplotype_variants_output_temp',num2str(h),');']);
                end
            end
            
            eval(['haplotype_analysis_output_final = haplotype_analysis_output',num2str(num_iterate),';']);
            writetable(haplotype_analysis_output_final,['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'_haplotype.txt'],'Delimiter','\t');
        end
        
        if perform_variant_analysis == 1 && end_number == starting_number && b==length(PAM_sequence_list) && output_to_excel == 0 && exist('merge_unique_variant_guides_full','var')==1 || perform_variant_analysis == 1 && individual_output == 1 && output_to_excel == 0 && b==length(PAM_sequence_list) && exist('merge_unique_variant_guides_full','var')==1 || perform_wgs_analysis == 1 && end_number == starting_number && b==length(PAM_sequence_list) && output_to_excel == 0 && exist('merge_unique_variant_guides_full','var')==1 || perform_wgs_analysis == 1 && individual_output == 1 && output_to_excel == 0 && b==length(PAM_sequence_list) && exist('merge_unique_variant_guides_full','var')==1;
            array_cell_variant = {'ID' 'PAM' 'sgRNA_Sequence' 'sgRNA_PAM' 'Full_Sequence' 'Strand' 'Chromosome' 'Classification' 'N_W_Alignment_Score' 'No_Variant_ID' 'No_Variant_sgRNA_Sequence' 'No_Variant_PAM_Only' 'No_Variant_Full_Sequence' 'No_Variant_Strand' 'No_Variant_Chromosome' 'Guide_Start_Coordinate' 'Guide_End_Coordinate' 'DSB_Coordinate' 'Guide_Start_Index' 'Guide_End_Index' 'DSB_Index' 'Filename'};
            merge_variant_output = table(merge_unique_variant_guides_full.guide_ID_number,merge_unique_variant_guides_full.pam,merge_unique_variant_guides_full.guides_only,merge_unique_variant_guides_full.pam_only,merge_unique_variant_guides_full.guides,merge_unique_variant_guides_full.strand,merge_unique_variant_guides_full.chr_id,merge_unique_variant_guides_full.pam_create,merge_unique_variant_guides_full.n_mm_wt,merge_unique_variant_guides_full.guide_ID_number_wt,merge_unique_variant_guides_full.guides_only_wt,merge_unique_variant_guides_full.pam_only_wt,merge_unique_variant_guides_full.guides_wt,merge_unique_variant_guides_full.strand_wt,merge_unique_variant_guides_full.chr_id_wt,merge_unique_variant_guides_full.guide_coord_start_wt,merge_unique_variant_guides_full.guide_coord_end_wt,merge_unique_variant_guides_full.dsb_coord_wt,merge_unique_variant_guides_full.guide_index_start_wt,merge_unique_variant_guides_full.guide_index_end_wt,merge_unique_variant_guides_full.dsb_index_wt,merge_unique_variant_guides_full.filename,'VariableNames',array_cell_variant);
            
            num_iterate = num_associated_variants_batch(a,1);
            counter = 1;
            for h = 1:num_iterate
                array_cell_variants(1,(counter:(counter+8))) = {['RS_ID_',num2str(h),''] ['Ref_SNP_',num2str(h),''] ['Alt_SNP_',num2str(h),''] ['Variant_Class_',num2str(h),''] ['Variant_Chromosome_',num2str(h),''] ['Variant_Start_Coordinate_',num2str(h),''] ['Variant_End_Coordinate_',num2str(h),''] ['Variant_Start_Index_',num2str(h),''] ['Variant_End_Index_',num2str(h),'']};
                counter = counter + 9;
            end
            counter = 1;
            for h = 1:num_iterate;
                eval(['variants_output_temp',num2str(h),' = table(merge_unique_variant_guides_full',num2str(h),'.rs_id_var,merge_unique_variant_guides_full',num2str(h),'.ref_snp_var,merge_unique_variant_guides_full',num2str(h),'.alt_snp_var,merge_unique_variant_guides_full',num2str(h),'.identifier_var,merge_unique_variant_guides_full',num2str(h),'.chr_id_var,merge_unique_variant_guides_full',num2str(h),'.bpstart_var,merge_unique_variant_guides_full',num2str(h),'.bpend_var,merge_unique_variant_guides_full',num2str(h),'.variant_index_start_var,merge_unique_variant_guides_full',num2str(h),'.variant_index_end_var,''VariableNames'',array_cell_variants(1,(counter:counter+8)));']);
                counter = counter + 9;
            end
            
            for h = 1:num_iterate;
                if h == 1;
                    eval(['variant_analysis_output',num2str(h),' = horzcat(merge_variant_output,variants_output_temp',num2str(h),');']);
                elseif h == 2
                    eval(['variant_analysis_output',num2str(h),' = horzcat(variant_analysis_output',num2str(h-1),',variants_output_temp',num2str(h),');']);
                elseif h>2
                    eval(['variant_analysis_output',num2str(h),' = horzcat(variant_analysis_output',num2str(h-1),',variants_output_temp',num2str(h),');']);
                end
            end
            
            eval(['variant_analysis_output_final = variant_analysis_output',num2str(num_iterate),';']);
            writetable(variant_analysis_output_final,['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'_variants.txt'],'Delimiter','\t');
        end
        
        if perform_haplotype_analysis == 1 && end_number == starting_number && output_to_excel == 1 && b==length(PAM_sequence_list) || perform_haplotype_analysis == 1 && individual_output == 1 && output_to_excel == 1 && b==length(PAM_sequence_list);
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'ID'}, 'Haplotype Analysis sgRNA', char(alphabet1(1)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_unique_haplotype_guides_full.guide_ID_number, 'Haplotype Analysis sgRNA', char(alphabet2(1)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'PAM'}, 'Haplotype Analysis sgRNA', char(alphabet1(2)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_unique_haplotype_guides_full.pam, 'Haplotype Analysis sgRNA', char(alphabet2(2)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'sgRNA Sequence'}, 'Haplotype Analysis sgRNA', char(alphabet1(3)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_unique_haplotype_guides_full.guides_only, 'Haplotype Analysis sgRNA', char(alphabet2(3)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'PAM'}, 'Haplotype Analysis sgRNA', char(alphabet1(4)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_unique_haplotype_guides_full.pam_only, 'Haplotype Analysis sgRNA', char(alphabet2(4)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'sgRNA+PAM'}, 'Haplotype Analysis sgRNA', char(alphabet1(5)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_unique_haplotype_guides_full.guides, 'Haplotype Analysis sgRNA', char(alphabet2(5)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'Strand'}, 'Haplotype Analysis sgRNA', char(alphabet1(6)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_unique_haplotype_guides_full.strand, 'Haplotype Analysis sgRNA', char(alphabet2(6)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'Chromosome'}, 'Haplotype Analysis sgRNA', char(alphabet1(7)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_unique_haplotype_guides_full.chr_id, 'Haplotype Analysis sgRNA', char(alphabet2(7)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'Count'}, 'Haplotype Analysis sgRNA', char(alphabet1(8)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_unique_haplotype_guides_full.count, 'Haplotype Analysis sgRNA', char(alphabet2(8)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'Frequency'}, 'Haplotype Analysis sgRNA', char(alphabet1(9)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_unique_haplotype_guides_full.frequency, 'Haplotype Analysis sgRNA', char(alphabet2(9)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'sgRNA Classification'}, 'Haplotype Analysis sgRNA', char(alphabet1(10)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_unique_haplotype_guides_full.pam_create, 'Haplotype Analysis sgRNA', char(alphabet2(10)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'N-W Alignment Score'}, 'Haplotype Analysis sgRNA', char(alphabet1(11)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_unique_haplotype_guides_full.n_mm_wt, 'Haplotype Analysis sgRNA', char(alphabet2(11)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'ID'}, 'Haplotype Analysis sgRNA', char(alphabet1(12)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_unique_haplotype_guides_full.guide_ID_number_wt, 'Haplotype Analysis sgRNA', char(alphabet2(12)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'Sequence of Associated No Variant sgRNA'}, 'Haplotype Analysis sgRNA', char(alphabet1(13)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_unique_haplotype_guides_full.guides_only_wt, 'Haplotype Analysis sgRNA', char(alphabet2(13)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'PAM'}, 'Haplotype Analysis sgRNA', char(alphabet1(14)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_unique_haplotype_guides_full.pam_only_wt, 'Haplotype Analysis sgRNA', char(alphabet2(14)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'sgRNA+PAM'}, 'Haplotype Analysis sgRNA', char(alphabet1(15)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_unique_haplotype_guides_full.guides_wt, 'Haplotype Analysis sgRNA', char(alphabet2(15)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'Strand'}, 'Haplotype Analysis sgRNA', char(alphabet1(16)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_unique_haplotype_guides_full.strand_wt, 'Haplotype Analysis sgRNA', char(alphabet2(16)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'Chromosome'}, 'Haplotype Analysis sgRNA', char(alphabet1(17)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_unique_haplotype_guides_full.chr_id_wt, 'Haplotype Analysis sgRNA', char(alphabet2(17)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], temp_coord_string_start, 'Haplotype Analysis sgRNA', char(alphabet1(18)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_unique_haplotype_guides_full.guide_coord_start_wt, 'Haplotype Analysis sgRNA', char(alphabet2(18)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], temp_coord_string_end, 'Haplotype Analysis sgRNA', char(alphabet1(19)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_unique_haplotype_guides_full.guide_coord_end_wt, 'Haplotype Analysis sgRNA', char(alphabet2(19)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], temp_coord_string_dsb, 'Haplotype Analysis sgRNA', char(alphabet1(20)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_unique_haplotype_guides_full.dsb_coord_wt, 'Haplotype Analysis sgRNA', char(alphabet2(20)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'Guide Start (Index)'}, 'Haplotype Analysis sgRNA', char(alphabet1(21)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_unique_haplotype_guides_full.guide_index_start_wt, 'Haplotype Analysis sgRNA', char(alphabet2(21)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'Guide End (Index)'}, 'Haplotype Analysis sgRNA', char(alphabet1(22)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_unique_haplotype_guides_full.guide_index_end_wt, 'Haplotype Analysis sgRNA', char(alphabet2(22)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'Double Strand Break (Index)'}, 'Haplotype Analysis sgRNA', char(alphabet1(23)));
            xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], merge_unique_haplotype_guides_full.dsb_index_wt, 'Haplotype Analysis sgRNA', char(alphabet2(23)));
            
            num_iterate = num_associated_variants_batch(a,1);
            
            for h = 1:num_iterate;
                xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'RS ID'}, 'Haplotype Analysis sgRNA', char(alphabet1(24+9*(h-1))));
                xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], eval(['merge_unique_haplotype_guides_full',num2str(h),'.rs_id_var']), 'Haplotype Analysis sgRNA', char(alphabet2(24+9*(h-1))));
                xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'Ref SNP'}, 'Haplotype Analysis sgRNA', char(alphabet1(25+9*(h-1))));
                xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], eval(['merge_unique_haplotype_guides_full',num2str(h),'.ref_snp_var']), 'Haplotype Analysis sgRNA', char(alphabet2(25+9*(h-1))));
                xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'Alt SNP'}, 'Haplotype Analysis sgRNA', char(alphabet1(26+9*(h-1))));
                xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], eval(['merge_unique_haplotype_guides_full',num2str(h),'.alt_snp_var']), 'Haplotype Analysis sgRNA', char(alphabet2(26+9*(h-1))));
                xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'Variant Class'}, 'Haplotype Analysis sgRNA', char(alphabet1(27+9*(h-1))));
                xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], eval(['merge_unique_haplotype_guides_full',num2str(h),'.identifier_var']), 'Haplotype Analysis sgRNA', char(alphabet2(27+9*(h-1))));
                xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'Variant Chr'}, 'Haplotype Analysis sgRNA', char(alphabet1(28+9*(h-1))));
                xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], eval(['merge_unique_haplotype_guides_full',num2str(h),'.chr_id_var']), 'Haplotype Analysis sgRNA', char(alphabet2(28+9*(h-1))));
                xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], temp_coord_string_start_variant, 'Haplotype Analysis sgRNA', char(alphabet1(29+9*(h-1))));
                xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], eval(['merge_unique_haplotype_guides_full',num2str(h),'.bpstart_var']), 'Haplotype Analysis sgRNA', char(alphabet2(29+9*(h-1))));
                xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], temp_coord_string_end_variant, 'Haplotype Analysis sgRNA', char(alphabet1(30+9*(h-1))));
                xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], eval(['merge_unique_haplotype_guides_full',num2str(h),'.bpend_var']), 'Haplotype Analysis sgRNA', char(alphabet2(30+9*(h-1))));
                xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'Variant Start (Index)'}, 'Haplotype Analysis sgRNA', char(alphabet1(31+9*(h-1))));
                xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], eval(['merge_unique_haplotype_guides_full',num2str(h),'.variant_index_start_var']), 'Haplotype Analysis sgRNA', char(alphabet2(31+9*(h-1))));
                xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], {'Variant End (Index)'}, 'Haplotype Analysis sgRNA', char(alphabet1(32+9*(h-1))));
                xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'.xlsx'], eval(['merge_unique_haplotype_guides_full',num2str(h),'.variant_index_end_var']), 'Haplotype Analysis sgRNA', char(alphabet2(32+9*(h-1))));
            end
        end
        
        if perform_plot_analysis == 1 && end_number == starting_number && b==length(PAM_sequence_list) || perform_plot_analysis == 1 && individual_output == 1 && b==length(PAM_sequence_list)
            
            p1 = 50;
            p2 = 100;
            p3 = 1200;
            p4 = 450;
            
            figure('position', [p1, p2, p3, p4])
            if perform_haplotype_analysis == 1
                haplo_guide_freq_index = find(merge_unique_haplotype_guides_full.frequency>haplo_guide_freq_cutoff);
                
            end
            
            pts = seq_bpstart:1:seq_bpend;
            [f,xi] = ksdensity(merge_original_guides_full.dsb_coord,pts,'width',ks_bandwidth);
            A(1:length(merge_original_guides_full.strand)) = (max(f)+min(f))/2;
            B(1:length(merge_original_guides_full.strand)) = (max(f)+min(f))/2;
            A = transpose(A);
            B = transpose(B);
            haplotype_change1 = (max(f)+min(f))*0.05;
            haplotype_change2 = (max(f)+min(f))*0.1;
            
            if perform_haplotype_analysis == 1
                E(1:length(haplo_guide_freq_index),1) = ((max(f)+min(f))/2)+haplotype_change2;
            end
            if perform_variant_analysis == 1 || perform_wgs_analysis == 1
                E(1:length(merge_unique_variant_guides_full.dsb_coord_wt),1) = ((max(f)+min(f))/2)+haplotype_change2;
            end
            G = [(((max(f)+min(f))/2)+haplotype_change1) (((max(f)+min(f))/2)+haplotype_change1)];
            
            markersize = 15;
            [ax, h(1), h(2)] = plotyy(merge_original_guides_full.dsb_coord,A,xi,f,'plot');
            set(h(1),'linestyle','none','marker','.','MarkerSize',markersize,'Color','k');
            set(h(2),'Color','k')
            hold on
            if perform_haplotype_analysis == 1
                if length(haplo_guide_freq_index)>0
                    h(6) = plot(merge_unique_haplotype_guides_full.dsb_coord_wt(haplo_guide_freq_index),E,'Parent', ax(1));
                    set(h(6),'linestyle','none','marker','.','MarkerSize',markersize,'Color','b');
                    
                end
                for j = 1:length(full_list_haplo_variants.bpstart)
                    h(8) = plot([full_list_haplo_variants.bpstart(j) full_list_haplo_variants.bpend(j)],G,'.','Parent', ax(1));
                    set(h(8),'linestyle','none','marker','.','MarkerSize',markersize,'Color','m');
                end
            end
            if perform_variant_analysis == 1 || perform_wgs_analysis == 1
                h(6) = plot(merge_unique_variant_guides_full.dsb_coord_wt,E,'Parent', ax(1));
                set(h(6),'linestyle','none','marker','.','MarkerSize',markersize,'Color','b');
                for j = 1:length(test_variants.bpstart)
                    h(8) = plot([test_variants.bpstart(j) test_variants.bpend(j)],G,'.','Parent', ax(1));
                    set(h(8),'linestyle','none','marker','.','MarkerSize',markersize,'Color','m');
                end
            end
            
            hold off
            set(gca,'ytick',[])
            set(ax,'ytick',[])
            set(gca,'ycolor',[1 1 1])
            set(ax,'ycolor',[1 1 1])
            set(ax,'YLim', [min(f) max(f)]);
            set(gca,'YLim', [min(f) max(f)]);
            xlabel(['Genomic Position (',coordinate_sys,')'])
            set(ax, 'XLim', [seq_bpstart seq_bpend]);
            set(gca, 'Xlim', [seq_bpstart seq_bpend]);
            customtick_top = get(gca, 'XTick');
            for n = 1:length(customtick_top)
                customtick_bot(1,n) = find(customtick_top(1,n)==reference_sequence.coord);
            end
            set(gca, 'XTickLabel', cellstr(num2str(customtick_top(:))));
            set(ax, 'XTickLabel', cellstr(num2str(customtick_top(:))));
            A = cellstr(num2str(customtick_top(:)));
            B = cellstr(num2str(customtick_bot(:)));
            if perform_haplotype_analysis == 1 && length(haplo_guide_freq_index)>0;
                [legh,objh,outh,outm] = legend([h(2),h(1),h(8),h(6)],['DSB density (KS function, bandwidth ',num2str(ks_bandwidth),')'],'sgRNA DSB','Haplotype variants',['Haplotype-associated sgRNA DSB, >',num2str(haplo_guide_freq_cutoff),'% guide frequency']);
                set(legh,'Location','SouthOutside');
                set(legh,'color','w')
                set(objh,'linewidth',1);
            elseif perform_haplotype_analysis == 1 && length(haplo_guide_freq_index)==0;
                [legh,objh,outh,outm] = legend([h(2),h(1),h(8)],'KS Density Function','sgRNA DSB','Haplotype Variants');
                set(legh,'Location','SouthOutside');
                set(legh,'color','w')
                set(objh,'linewidth',1);
            end
            if perform_variant_analysis == 1 || perform_wgs_analysis == 1
                if perform_wgs_analysis == 0 && perform_variant_analysis == 1
                    [legh,objh,outh,outm] = legend([h(2),h(1),h(8),h(6)],['DSB density (KS function, bandwidth ',num2str(ks_bandwidth),')'],'sgRNA DSB','Custom list of variants','Variant-associated sgRNA DSB');
                elseif perform_wgs_analysis == 1 && perform_variant_analysis == 0
                    [legh,objh,outh,outm] = legend([h(2),h(1),h(8),h(6)],['DSB density (KS function, bandwidth ',num2str(ks_bandwidth),')'],'sgRNA DSB','WGS variants','WGS-associated sgRNA DSB');
                end
                set(legh,'Location','SouthOutside');
                set(legh,'color','w')
                set(objh,'linewidth',1);
            end
            
            if perform_haplotype_analysis == 0 && perform_variant_analysis == 0 && perform_wgs_analysis == 0
                [legh,objh,outh,outm] = legend([h(2),h(1)],['DSB density (KS function, bandwidth ',num2str(ks_bandwidth),')'],'sgRNA DSB');
                set(legh,'Location','SouthOutside');
                set(legh,'color','w')
                set(objh,'linewidth',1);
            end
            
            if save_figures == 1;
                h=gcf;
                set(h,'PaperPositionMode','auto');
                set(h,'PaperOrientation','landscape');
                set(h,'Position',[p1 p2 p3 p4]);
                print(gcf, '-dpdf', ['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'_Fig1.pdf']);
                print(gcf,'-dtiff','-r300',['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'_Fig1.tiff'])
            end
        end
        
        if b == length(PAM_sequence_list)
            dsb_index_temp = sort(merge_original_guides_full.dsb_index,1);
            sgRNA_gaps.gap_distance = diff(dsb_index_temp);
            max_gap = max(sgRNA_gaps.gap_distance);
            percentile_50(a,1) = prctile(sgRNA_gaps.gap_distance,50);
            percentile_90(a,1) = prctile(sgRNA_gaps.gap_distance,90);
            if length(max_gap)>0
                max_gap_list(a,1) = max_gap;
            else
                max_gap_list(a,1) = 0;
            end
            
            if individual_output == 1 && output_to_excel == 1
                xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'_gaps.xlsx'],{'Filename'},'Gaps','A1');
                if length(ref_seq_info.identifier)>=1;
                    xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'_gaps.xlsx'],ref_seq_info.identifier(a,1),'Gaps','A2');
                end
                xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'_gaps.xlsx'],{'50th Percentile Gap Distance'},'Gaps','B1');
                xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'_gaps.xlsx'],percentile_50(a,1),'Gaps','B2');
                xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'_gaps.xlsx'],{'90th Percentile Gap Distance'},'Gaps','C1');
                xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'_gaps.xlsx'],percentile_90(a,1),'Gaps','C2');
                xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'_gaps.xlsx'],{'Maximum Gap Distance'},'Gaps','D1');
                xlswrite(['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'_gaps.xlsx'],max_gap_list(a,1),'Gaps','D2');
            end
            if individual_output == 1 && output_to_excel == 0
                if length(ref_seq_info.identifier) >= 1
                    array_cell_no_variant_gaps = {'Filename' 'Fiftieth_Percentile_Gap_Distance' 'Ninetieth_Percentile_Gap_Distance' 'Maximum_Gap_Distance'};
                    no_variant_output_gaps = table(ref_seq_info.identifier(a,1),percentile_50(a,1),percentile_90(a,1),max_gap_list(a,1),'VariableNames',array_cell_no_variant_gaps);
                    writetable(no_variant_output_gaps,['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'_no_variant_gaps.txt'],'Delimiter','\t')
                elseif length(ref_seq_info.identifier) < 1
                    array_cell_no_variant_gaps = {'Fiftieth_Percentile_Gap_Distance' 'Ninetieth_Percentile_Gap_Distance' 'Maximum_Gap_Distance'};
                    no_variant_output_gaps = table(percentile_50(a,1),percentile_90(a,1),max_gap_list(a,1),'VariableNames',array_cell_no_variant_gaps);
                    writetable(no_variant_output_gaps,['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'_no_variant_gaps.txt'],'Delimiter','\t')
                end
            end
        end
        
        if perform_plot_analysis == 1 && end_number > starting_number && perform_batch_output == 1 && exist('sgRNA_gaps','var')==1 && b == length(PAM_sequence_list)
            sgRNA_gaps_full.gap_distances = [sgRNA_gaps_full.gap_distances;sgRNA_gaps.gap_distance];
            sgRNA_gaps_full_filename_temp(1:length(sgRNA_gaps.gap_distance)) = cellstr(temp_filename);
            sgRNA_gaps_full.filenames = [sgRNA_gaps_full.filenames;transpose(sgRNA_gaps_full_filename_temp)];
        end
        
        if perform_plot_analysis == 1 && end_number == starting_number && exist('sgRNA_gaps','var')==1 && b == length(PAM_sequence_list) || perform_plot_analysis == 1 && individual_output == 1 && exist('sgRNA_gaps','var')==1 && b == length(PAM_sequence_list)
            [histo_counts1 histo_centers1] = hist(sgRNA_gaps.gap_distance,length(sgRNA_gaps.gap_distance));
            
            histo_counts_index = find(histo_counts1>0);
            histo_counts = transpose(histo_counts1(histo_counts_index));
            histo_centers = transpose(round(histo_centers1(histo_counts_index),0));
            histo_fraction_total = sum(histo_counts);
            histo_fraction = histo_counts/histo_fraction_total;
            histo_fraction_percent = 100*histo_fraction;
            
            figure;
            bar(histo_centers,histo_fraction_percent);
            y1=get(gca,'ylim');
            hold on
            h1 = plot([prctile(sgRNA_gaps.gap_distance,50) prctile(sgRNA_gaps.gap_distance,50)],y1,'m');
            h2 = plot([prctile(sgRNA_gaps.gap_distance,90) prctile(sgRNA_gaps.gap_distance,90)],y1,'r');
            hold off
            xlim([-2 max(sgRNA_gaps.gap_distance)])
            xlabel('Distance between adjacent genomic cleavages (bp)')
            ylabel('Percent (%)');
            [legh,objh,outh,outm] = legend([h1 h2],{['50^{th} percentile (',num2str(prctile(sgRNA_gaps.gap_distance,50)),' bp)'],['90^{th} percentile (',num2str(prctile(sgRNA_gaps.gap_distance,90)),' bp)']});
            set(legh,'Location','NorthOutside');
            set(legh,'color','w')
            
            if save_figures == 1;
                h=gcf;
                set(h,'PaperPositionMode','auto');
                set(h,'PaperOrientation','landscape');
                print(gcf, '-dpdf', ['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'_Fig2.pdf']);
                print(gcf,'-dtiff','-r300',['OUTPUT_',num2str(a),'_',reference_seq_title,'_PAM=',PAM_seq_output,'_',coordinate_sys,'_Fig2.tiff']);
            end
        end
        
        if end_number>starting_number || length(PAM_sequence_list)>1
            clearvars -except percentile_50 percentile_90 batch_filename_gaps max_gap_list perform_wgs_analysis batch_unique_variant_guides_full merge_unique_variant_guides_full num_spots_batch num_spots_sum num_associated_variants_batch PAM_seq_output ks_bandwidth batch_multiple_match multiple_match_analysis a PAM_sequence_list sgRNA_len_list five_prime_PAM_list output_to_excel temp_coord_string_start_variant temp_coord_string_end_variant num_associated_variants num_spots ks_bandwidth batch_filename_figure save_figures haplo_guide_freq_cutoff alphabet1 alphabet2 perform_plot_analysis plot_dsb_position perform_variant_analysis five_prime_PAM twenty_extra_bases sgRNA_gaps_full batch_filename batch_original_guides_full batch_unique_haplotype_guides_full individual_output perform_haplotype_analysis max_number_guide_ID_digits algorithm_type algorithm_type2 perform_batch_output coordinate_sys starting_number window_size end_number reference_seq1 ref_seq_info vcf_filenames reference_seq_title_file ref_sequence_info Header num_seq_upload guide_ID_num_temp vcf_haplo_filenames temp_coord_string_start temp_coord_string_end temp_coord_string_dsb no_variant_output_num_col haplotype_outout_num_col output_filenames_list batch_unique_haplotype_guides_full1 batch_unique_haplotype_guides_full2 batch_unique_haplotype_guides_full3 batch_unique_haplotype_guides_full4 batch_unique_haplotype_guides_full5 batch_unique_haplotype_guides_full6 batch_unique_haplotype_guides_full7 batch_unique_haplotype_guides_full8 merge_multiple_match merge_original_guides_full merge_unique_haplotype_guides_full merge_unique_haplotype_guides_full1 merge_unique_haplotype_guides_full2 merge_unique_haplotype_guides_full3 merge_unique_haplotype_guides_full4 merge_unique_haplotype_guides_full5 merge_unique_haplotype_guides_full6 merge_unique_haplotype_guides_full7 merge_unique_haplotype_guides_full8 merge_unique_variant_guides_full1 merge_unique_variant_guides_full2 merge_unique_variant_guides_full3 merge_unique_variant_guides_full4 merge_unique_variant_guides_full5 merge_unique_variant_guides_full6 merge_unique_variant_guides_full7 merge_unique_variant_guides_full8 batch_unique_variant_guides_full1 batch_unique_variant_guides_full2 batch_unique_variant_guides_full3 batch_unique_variant_guides_full4 batch_unique_variant_guides_full5 batch_unique_variant_guides_full6 batch_unique_variant_guides_full7 batch_unique_variant_guides_full8
        end
    end
    if perform_haplotype_analysis == 1 && exist('num_associated_variants','var')==1 || perform_variant_analysis == 1 && exist('num_associated_variants','var')==1 || perform_wgs_analysis == 1 && exist('num_associated_variants','var')==1
        num_associated_variants_max(a,1) = max(num_associated_variants);
    end
    clearvars num_spots num_spots_sum num_associated_variants merge_multiple_match merge_original_guides_full merge_unique_haplotype_guides_full merge_unique_haplotype_guides_full1 merge_unique_haplotype_guides_full2 merge_unique_haplotype_guides_full3 merge_unique_haplotype_guides_full4 merge_unique_haplotype_guides_full5 merge_unique_haplotype_guides_full6 merge_unique_haplotype_guides_full7 merge_unique_haplotype_guides_full8 merge_unique_variant_guides_full merge_unique_variant_guides_full1 merge_unique_variant_guides_full2 merge_unique_variant_guides_full3 merge_unique_variant_guides_full4 merge_unique_variant_guides_full5 merge_unique_variant_guides_full6 merge_unique_variant_guides_full7 merge_unique_variant_guides_full8
end

%% Batch Analysis

if end_number > starting_number && perform_batch_output == 1 && perform_plot_analysis == 1
    [histo_counts1 histo_centers1] = hist(sgRNA_gaps_full.gap_distances,length(sgRNA_gaps_full.gap_distances));
    
    histo_counts_index = find(histo_counts1>0);
    histo_counts = transpose(histo_counts1(histo_counts_index));
    histo_centers = transpose(round(histo_centers1(histo_counts_index),0));
    histo_fraction_total = sum(histo_counts);
    histo_fraction = histo_counts/histo_fraction_total;
    histo_fraction_percent = 100*histo_fraction;
    
    figure;
    bar(histo_centers,histo_fraction_percent);
    y1=get(gca,'ylim');
    hold on
    h1 = plot([prctile(sgRNA_gaps_full.gap_distances,50) prctile(sgRNA_gaps_full.gap_distances,50)],y1,'m');
    h2 = plot([prctile(sgRNA_gaps_full.gap_distances,90) prctile(sgRNA_gaps_full.gap_distances,90)],y1,'r');
    hold off
    xlim([-2 max(sgRNA_gaps_full.gap_distances)])
    
    xlabel('Distance between adjacent genomic cleavages (bp)')
    ylabel('Percent (%)');
    [legh,objh,outh,outm] = legend([h1 h2],{['50^{th} percentile (',num2str(prctile(sgRNA_gaps_full.gap_distances,50)),' bp)'],['90^{th} percentile (',num2str(prctile(sgRNA_gaps_full.gap_distances,90)),' bp)']});
    set(legh,'Location','NorthOutside');
    set(legh,'color','w')
    
    if save_figures == 1;
        h=gcf;
        set(h,'PaperPositionMode','auto');
        set(h,'PaperOrientation','landscape');
        print(gcf, '-dpdf', batch_filename_figure)
    end
end

if end_number > starting_number && perform_batch_output == 1 && output_to_excel == 0
    guide_ID_num_temp = 1;
    for p = 1:length(batch_original_guides_full.guides)
        zeros_to_add = [];
        num_zeros_add = length(num2str(guide_ID_num_temp));
        if num_zeros_add < max_number_guide_ID_digits
            for j = 1:(max_number_guide_ID_digits-num_zeros_add)
                zeros_to_add = strcat(zeros_to_add,'0');
            end
            batch_original_guides_full.guide_ID_number{p,1} = strcat('sgRNA_',zeros_to_add,num2str(guide_ID_num_temp));
        elseif num_zeros_add == max_number_guide_ID_digits
            batch_original_guides_full.guide_ID_number{p,1} = ['sgRNA_',num2str(guide_ID_num_temp),''];
        end
        guide_ID_num_temp = guide_ID_num_temp + 1;
    end
    array_cell_no_variant = {'ID' 'PAM' 'sgRNA_Sequence' 'PAM_Only' 'sgRNA_PAM' 'Strand' 'Chromosome' 'Guide_Start_Coordinate' 'Guide_End_Coordinate' 'DSB_Coordinate' 'Guide_Start_Index' 'Guide_End_Index' 'Double_Strand_Break_Index' 'Filename'};
    batch_no_variant_output = table(batch_original_guides_full.guide_ID_number,batch_original_guides_full.pam,batch_original_guides_full.guides_only,batch_original_guides_full.pam_only,batch_original_guides_full.guides,batch_original_guides_full.strand,batch_original_guides_full.chr_id,batch_original_guides_full.guide_coord_start,batch_original_guides_full.guide_coord_end,batch_original_guides_full.dsb_coord,batch_original_guides_full.guide_index_start,batch_original_guides_full.guide_index_end,batch_original_guides_full.dsb_index,batch_original_guides_full.filename,'VariableNames',array_cell_no_variant);
    writetable(batch_no_variant_output,['OUTPUT_',num2str(a),'_PAM=',PAM_seq_output,'_',coordinate_sys,'_no_variant_batch.txt'],'Delimiter','\t')
    
    if length(ref_seq_info.identifier) >= 1
        array_cell_no_variant_gaps = {'Filename' 'Fiftieth_Percentile_Gap_Distance' 'Ninetieth_Percentile_Gap_Distance' 'Maximum_Gap_Distance'};
        batch_no_variant_output_gaps = table(ref_seq_info.identifier,percentile_50,percentile_90,max_gap_list,'VariableNames',array_cell_no_variant_gaps);
        writetable(batch_no_variant_output_gaps,['OUTPUT_',num2str(a),'_PAM=',PAM_seq_output,'_',coordinate_sys,'_no_variant_batch_gaps.txt'],'Delimiter','\t')
    elseif length(ref_seq_info.identifier) < 1
        array_cell_no_variant_gaps = {'Fiftieth_Percentile_Gap_Distance' 'Ninetieth_Percentile_Gap_Distance' 'Maximum_Gap_Distance'};
        batch_no_variant_output_gaps = table(percentile_50,percentile_90,max_gap_list,'VariableNames',array_cell_no_variant_gaps);
        writetable(batch_no_variant_output_gaps,['OUTPUT_',num2str(a),'_PAM=',PAM_seq_output,'_',coordinate_sys,'_no_variant_batch_gaps.txt'],'Delimiter','\t')
    end
end



if multiple_match_analysis == 1 && perform_batch_output == 1 && exist('batch_multiple_match','var')==1 && end_number>starting_number && output_to_excel == 0
    array_cell_multiple_match = {'ID' 'PAM' 'sgRNA_Sequence' 'PAM_Only' 'sgRNA_PAM' 'Strand' 'Chromosome' 'Guide_Start_Coordinate' 'Guide_End_Coordinate' 'DSB_Coordinate' 'Guide_Start_Index' 'Guide_End_Index' 'Double_Strand_Break_Index' 'Filename'};
    batch_multiple_match_output = table(batch_multiple_match.guide_ID_number,batch_multiple_match.pam,batch_multiple_match.guides_only,batch_multiple_match.pam_only,batch_multiple_match.guides,batch_multiple_match.strand,batch_multiple_match.chr_id,batch_multiple_match.guide_coord_start,batch_multiple_match.guide_coord_end,batch_multiple_match.dsb_coord,batch_multiple_match.guide_index_start,batch_multiple_match.guide_index_end,batch_multiple_match.dsb_index,batch_multiple_match.filename,'VariableNames',array_cell_multiple_match);
    writetable(batch_multiple_match_output,['OUTPUT_',num2str(a),'_PAM=',PAM_seq_output,'_',coordinate_sys,'_multiple_match_batch.txt'],'Delimiter','\t')
end

if perform_haplotype_analysis == 1 && end_number > starting_number && output_to_excel == 0 && perform_batch_output == 1
    guide_ID_num_temp = length(batch_original_guides_full.guides)+1;
    for p = 1:length(batch_unique_haplotype_guides_full.guides)
        zeros_to_add = [];
        num_zeros_add = length(num2str(guide_ID_num_temp));
        if num_zeros_add < max_number_guide_ID_digits
            for j = 1:(max_number_guide_ID_digits-num_zeros_add)
                zeros_to_add = strcat(zeros_to_add,'0');
            end
            batch_unique_haplotype_guides_full.guide_ID_number{p,1} = strcat('sgRNA_',zeros_to_add,num2str(guide_ID_num_temp));
        elseif num_zeros_add == max_number_guide_ID_digits
            batch_unique_haplotype_guides_full.guide_ID_number{p,1} = ['sgRNA_',num2str(guide_ID_num_temp),''];
        end
        guide_ID_num_temp = guide_ID_num_temp + 1;
    end
    
    for p = 1:length(batch_unique_haplotype_guides_full.guides)
        renumb_temp = strmatch(batch_unique_haplotype_guides_full.guides_wt(p,1),batch_original_guides_full.guides);
        if renumb_temp > 0
            batch_unique_haplotype_guides_full.guide_ID_number_wt(p,1) =  batch_original_guides_full.guide_ID_number(renumb_temp);
        else
            batch_unique_haplotype_guides_full.guide_ID_number_wt(p,1) = {'N/A'};
        end
    end
    
    array_cell_haplotype = {'ID' 'PAM' 'sgRNA_Sequence' 'sgRNA_PAM' 'Full_Sequence' 'Strand' 'Chromosome' 'Count' 'Frequency' 'Classification' 'N_W_Alignment_Score' 'No_Variant_ID' 'No_Variant_sgRNA_Sequence' 'No_Variant_PAM_Only' 'No_Variant_Full_Sequence' 'No_Variant_Strand' 'No_Variant_Chromosome' 'Guide_Start_Coordinate' 'Guide_End_Coordinate' 'DSB_Coordinate' 'Guide_Start_Index' 'Guide_End_Index' 'DSB_Index' 'Filename'};
    batch_haplotype_output = table(batch_unique_haplotype_guides_full.guide_ID_number,batch_unique_haplotype_guides_full.pam,batch_unique_haplotype_guides_full.guides_only,batch_unique_haplotype_guides_full.pam_only,batch_unique_haplotype_guides_full.guides,batch_unique_haplotype_guides_full.strand,batch_unique_haplotype_guides_full.chr_id,batch_unique_haplotype_guides_full.count,batch_unique_haplotype_guides_full.frequency,batch_unique_haplotype_guides_full.pam_create,batch_unique_haplotype_guides_full.n_mm_wt,batch_unique_haplotype_guides_full.guide_ID_number_wt,batch_unique_haplotype_guides_full.guides_only_wt,batch_unique_haplotype_guides_full.pam_only_wt,batch_unique_haplotype_guides_full.guides_wt,batch_unique_haplotype_guides_full.strand_wt,batch_unique_haplotype_guides_full.chr_id_wt,batch_unique_haplotype_guides_full.guide_coord_start_wt,batch_unique_haplotype_guides_full.guide_coord_end_wt,batch_unique_haplotype_guides_full.dsb_coord_wt,batch_unique_haplotype_guides_full.guide_index_start_wt,batch_unique_haplotype_guides_full.guide_index_end_wt,batch_unique_haplotype_guides_full.dsb_index_wt,batch_unique_haplotype_guides_full.filename,'VariableNames',array_cell_haplotype);
    
    num_iterate = max(num_associated_variants_batch(a,1));
    counter = 1;
    for h = 1:num_iterate
        array_cell_haplotype_variants(1,(counter:(counter+8))) = {['RS_ID_',num2str(h),''] ['Ref_SNP_',num2str(h),''] ['Alt_SNP_',num2str(h),''] ['Variant_Class_',num2str(h),''] ['Variant_Chromosome_',num2str(h),''] ['Variant_Start_Coordinate_',num2str(h),''] ['Variant_End_Coordinate_',num2str(h),''] ['Variant_Start_Index_',num2str(h),''] ['Variant_End_Index_',num2str(h),'']};
        counter = counter + 9;
    end
    counter = 1;
    for h = 1:num_iterate;
        eval(['batch_haplotype_variants_output_temp',num2str(h),' = table(batch_unique_haplotype_guides_full',num2str(h),'.rs_id_var,batch_unique_haplotype_guides_full',num2str(h),'.ref_snp_var,batch_unique_haplotype_guides_full',num2str(h),'.alt_snp_var,batch_unique_haplotype_guides_full',num2str(h),'.identifier_var,batch_unique_haplotype_guides_full',num2str(h),'.chr_id_var,batch_unique_haplotype_guides_full',num2str(h),'.bpstart_var,batch_unique_haplotype_guides_full',num2str(h),'.bpend_var,batch_unique_haplotype_guides_full',num2str(h),'.variant_index_start_var,batch_unique_haplotype_guides_full',num2str(h),'.variant_index_end_var,''VariableNames'',array_cell_haplotype_variants(1,(counter:counter+8)));']);
        counter = counter + 9;
    end
    for h = 1:num_iterate;
        if h == 1;
            eval(['batch_haplotype_analysis_output',num2str(h),' = horzcat(batch_haplotype_output,batch_haplotype_variants_output_temp',num2str(h),');']);
        elseif h == 2
            eval(['batch_haplotype_analysis_output',num2str(h),' = horzcat(batch_haplotype_analysis_output',num2str(h-1),',batch_haplotype_variants_output_temp',num2str(h),');']);
        elseif h>2
            eval(['batch_haplotype_analysis_output',num2str(h),' = horzcat(batch_haplotype_analysis_output',num2str(h-1),',batch_haplotype_variants_output_temp',num2str(h),');']);
        end
    end
    
    eval(['batch_haplotype_analysis_output_final = batch_haplotype_analysis_output',num2str(num_iterate),';']);
    writetable(batch_haplotype_analysis_output_final,['OUTPUT_',num2str(a),'_PAM=',PAM_seq_output,'_',coordinate_sys,'_haplotype_batch.txt'],'Delimiter','\t');
end

if perform_variant_analysis == 1 && end_number > starting_number && output_to_excel == 0 && perform_batch_output == 1 && length(batch_unique_variant_guides_full.guides)>0 || perform_wgs_analysis == 1 && end_number > starting_number && output_to_excel == 0 && perform_batch_output == 1 && length(batch_unique_variant_guides_full.guides)>0
    guide_ID_num_temp = length(batch_original_guides_full.guides)+1;
    for p = 1:length(batch_unique_variant_guides_full.guides)
        zeros_to_add = [];
        num_zeros_add = length(num2str(guide_ID_num_temp));
        if num_zeros_add < max_number_guide_ID_digits
            for j = 1:(max_number_guide_ID_digits-num_zeros_add)
                zeros_to_add = strcat(zeros_to_add,'0');
            end
            batch_unique_variant_guides_full.guide_ID_number{p,1} = strcat('sgRNA_',zeros_to_add,num2str(guide_ID_num_temp));
        elseif num_zeros_add == max_number_guide_ID_digits
            batch_unique_variant_guides_full.guide_ID_number{p,1} = ['sgRNA_',num2str(guide_ID_num_temp),''];
        end
        guide_ID_num_temp = guide_ID_num_temp + 1;
    end
    
    for p = 1:length(batch_unique_variant_guides_full.guides)
        renumb_temp = strmatch(batch_unique_variant_guides_full.guides_wt(p,1),batch_original_guides_full.guides);
        if renumb_temp > 0
            batch_unique_variant_guides_full.guide_ID_number_wt(p,1) =  batch_original_guides_full.guide_ID_number(renumb_temp);
        else
            batch_unique_variant_guides_full.guide_ID_number_wt(p,1) = {'N/A'};
        end
    end
         
    array_cell_variant = {'ID' 'PAM' 'sgRNA_Sequence' 'sgRNA_PAM' 'Full_Sequence' 'Strand' 'Chromosome' 'Classification' 'N_W_Alignment_Score' 'No_Variant_ID' 'No_Variant_sgRNA_Sequence' 'No_Variant_PAM_Only' 'No_Variant_Full_Sequence' 'No_Variant_Strand' 'No_Variant_Chromosome' 'Guide_Start_Coordinate' 'Guide_End_Coordinate' 'DSB_Coordinate' 'Guide_Start_Index' 'Guide_End_Index' 'DSB_Index' 'Filename'};
    batch_variant_output = table(batch_unique_variant_guides_full.guide_ID_number,batch_unique_variant_guides_full.pam,batch_unique_variant_guides_full.guides_only,batch_unique_variant_guides_full.pam_only,batch_unique_variant_guides_full.guides,batch_unique_variant_guides_full.strand,batch_unique_variant_guides_full.chr_id,batch_unique_variant_guides_full.pam_create,batch_unique_variant_guides_full.n_mm_wt,batch_unique_variant_guides_full.guide_ID_number_wt,batch_unique_variant_guides_full.guides_only_wt,batch_unique_variant_guides_full.pam_only_wt,batch_unique_variant_guides_full.guides_wt,batch_unique_variant_guides_full.strand_wt,batch_unique_variant_guides_full.chr_id_wt,batch_unique_variant_guides_full.guide_coord_start_wt,batch_unique_variant_guides_full.guide_coord_end_wt,batch_unique_variant_guides_full.dsb_coord_wt,batch_unique_variant_guides_full.guide_index_start_wt,batch_unique_variant_guides_full.guide_index_end_wt,batch_unique_variant_guides_full.dsb_index_wt,batch_unique_variant_guides_full.filename,'VariableNames',array_cell_variant);
    
    if length(num_associated_variants_batch)<(end_number-starting_number+1)
        len_batch_temp = length(num_associated_variants_batch)+1;
        num_associated_variants_batch(len_batch_temp:end_number,1) = 0;
    end
    
    num_iterate = max(num_associated_variants_batch);
    counter = 1;
    for h = 1:num_iterate
        array_cell_variants(1,(counter:(counter+8))) = {['RS_ID_',num2str(h),''] ['Ref_SNP_',num2str(h),''] ['Alt_SNP_',num2str(h),''] ['Variant_Class_',num2str(h),''] ['Variant_Chromosome_',num2str(h),''] ['Variant_Start_Coordinate_',num2str(h),''] ['Variant_End_Coordinate_',num2str(h),''] ['Variant_Start_Index_',num2str(h),''] ['Variant_End_Index_',num2str(h),'']};
        counter = counter + 9;
    end
    counter = 1;
    for h = 1:num_iterate;
        eval(['batch_variants_output_temp',num2str(h),' = table(batch_unique_variant_guides_full',num2str(h),'.rs_id_var,batch_unique_variant_guides_full',num2str(h),'.ref_snp_var,batch_unique_variant_guides_full',num2str(h),'.alt_snp_var,batch_unique_variant_guides_full',num2str(h),'.identifier_var,batch_unique_variant_guides_full',num2str(h),'.chr_id_var,batch_unique_variant_guides_full',num2str(h),'.bpstart_var,batch_unique_variant_guides_full',num2str(h),'.bpend_var,batch_unique_variant_guides_full',num2str(h),'.variant_index_start_var,batch_unique_variant_guides_full',num2str(h),'.variant_index_end_var,''VariableNames'',array_cell_variants(1,(counter:counter+8)));']);
        counter = counter + 9;
    end
    for h = 1:num_iterate;
        if h == 1;
            eval(['batch_variant_analysis_output',num2str(h),' = horzcat(batch_variant_output,batch_variants_output_temp',num2str(h),');']);
        elseif h == 2
            eval(['batch_variant_analysis_output',num2str(h),' = horzcat(batch_variant_analysis_output',num2str(h-1),',batch_variants_output_temp',num2str(h),');']);
        elseif h>2
            eval(['batch_variant_analysis_output',num2str(h),' = horzcat(batch_variant_analysis_output',num2str(h-1),',batch_variants_output_temp',num2str(h),');']);
        end
    end
    eval(['batch_variant_analysis_output_final = batch_variant_analysis_output',num2str(num_iterate),';']);
    writetable(batch_variant_analysis_output_final,['OUTPUT_',num2str(a),'_PAM=',PAM_seq_output,'_',coordinate_sys,'_variants_batch.txt'],'Delimiter','\t');
end

if end_number > starting_number && perform_batch_output == 1 && output_to_excel == 1
    xlswrite(batch_filename, {'ID'}, 'No Variant Analysis', char(alphabet1(1)));
    xlswrite(batch_filename, batch_original_guides_full.guide_ID_number, 'No Variant Analysis', char(alphabet2(1)));
    xlswrite(batch_filename, {'PAM'}, 'No Variant Analysis', char(alphabet1(2)));
    xlswrite(batch_filename, batch_original_guides_full.pam, 'No Variant Analysis', char(alphabet2(2)));
    xlswrite(batch_filename, {'sgRNA Sequence'}, 'No Variant Analysis', char(alphabet1(3)));
    xlswrite(batch_filename, batch_original_guides_full.guides_only, 'No Variant Analysis', char(alphabet2(3)));
    xlswrite(batch_filename, {'PAM'}, 'No Variant Analysis', char(alphabet1(4)));
    xlswrite(batch_filename, batch_original_guides_full.pam_only, 'No Variant Analysis', char(alphabet2(4)));
    xlswrite(batch_filename, {'sgRNA+PAM'}, 'No Variant Analysis', char(alphabet1(5)));
    xlswrite(batch_filename, batch_original_guides_full.guides, 'No Variant Analysis', char(alphabet2(5)));
    xlswrite(batch_filename, {'Strand'}, 'No Variant Analysis', char(alphabet1(6)));
    xlswrite(batch_filename, batch_original_guides_full.strand, 'No Variant Analysis', char(alphabet2(6)));
    xlswrite(batch_filename, {'Chromosome'}, 'No Variant Analysis', char(alphabet1(7)));
    xlswrite(batch_filename, batch_original_guides_full.chr_id, 'No Variant Analysis', char(alphabet2(7)));
    xlswrite(batch_filename, temp_coord_string_start, 'No Variant Analysis', char(alphabet1(8)));
    xlswrite(batch_filename, batch_original_guides_full.guide_coord_start, 'No Variant Analysis', char(alphabet2(8)));
    xlswrite(batch_filename, temp_coord_string_end, 'No Variant Analysis', char(alphabet1(9)));
    xlswrite(batch_filename, batch_original_guides_full.guide_coord_end, 'No Variant Analysis', char(alphabet2(9)));
    xlswrite(batch_filename, temp_coord_string_dsb, 'No Variant Analysis', char(alphabet1(10)));
    xlswrite(batch_filename, batch_original_guides_full.dsb_coord, 'No Variant Analysis', char(alphabet2(10)));
    xlswrite(batch_filename, {'Guide Start (Index)'}, 'No Variant Analysis', char(alphabet1(11)));
    xlswrite(batch_filename, batch_original_guides_full.guide_index_start, 'No Variant Analysis', char(alphabet2(11)));
    xlswrite(batch_filename, {'Guide End (Index)'}, 'No Variant Analysis', char(alphabet1(12)));
    xlswrite(batch_filename, batch_original_guides_full.guide_index_end, 'No Variant Analysis', char(alphabet2(12)));
    xlswrite(batch_filename, {'Double Strand Break (Index)'}, 'No Variant Analysis', char(alphabet1(13)));
    xlswrite(batch_filename, batch_original_guides_full.dsb_index, 'No Variant Analysis', char(alphabet2(13)));
    if multiple_match_analysis == 1
        xlswrite(batch_filename, {'Predicted # of DSBs in Sequence'}, 'No Variant Analysis', char(alphabet1(14)));
        xlswrite(batch_filename, batch_original_guides_full.mult_match_count, 'No Variant Analysis', char(alphabet2(14)));
    end
    xlswrite(batch_filename, {'Filename'}, 'No Variant Analysis', char(alphabet1(15)));
    xlswrite(batch_filename, batch_original_guides_full.filename, 'No Variant Analysis', char(alphabet2(15)));
    
    xlswrite(batch_filename_gaps,{'Filename'},'Gaps','A1');
    if length(ref_seq_info.identifier)>=1;
        xlswrite(batch_filename_gaps,ref_seq_info.identifier,'Gaps','A2');
    end
    xlswrite(batch_filename_gaps,{'50th Percentile Gap Distance'},'Gaps','B1');
    xlswrite(batch_filename_gaps,percentile_50,'Gaps','B2');
    xlswrite(batch_filename_gaps,{'90th Percentile Gap Distance'},'Gaps','C1');
    xlswrite(batch_filename_gaps,percentile_90,'Gaps','C2');
    xlswrite(batch_filename_gaps,{'Maximum Gap Distance'},'Gaps','D1');
    xlswrite(batch_filename_gaps,max_gap_list,'Gaps','D2');
end

if end_number > starting_number && perform_batch_output == 1 && multiple_match_analysis == 1 && output_to_excel == 1 && length(batch_multiple_match.guides)>0
    xlswrite(batch_filename, {'ID'}, 'Multiple Match Analysis', char(alphabet1(1)));
    xlswrite(batch_filename, batch_multiple_match.guide_ID_number, 'Multiple Match Analysis', char(alphabet2(1)));
    xlswrite(batch_filename, {'PAM'}, 'Multiple Match Analysis', char(alphabet1(2)));
    xlswrite(batch_filename, batch_multiple_match.pam, 'Multiple Match Analysis', char(alphabet2(2)));
    xlswrite(batch_filename, {'sgRNA Sequence'}, 'Multiple Match Analysis', char(alphabet1(3)));
    xlswrite(batch_filename, batch_multiple_match.guides_only, 'Multiple Match Analysis', char(alphabet2(3)));
    xlswrite(batch_filename, {'PAM Only'}, 'Multiple Match Analysis', char(alphabet1(4)));
    xlswrite(batch_filename, batch_multiple_match.pam_only, 'Multiple Match Analysis', char(alphabet2(4)));
    xlswrite(batch_filename, {'sgRNA+PAM'}, 'Multiple Match Analysis', char(alphabet1(5)));
    xlswrite(batch_filename, batch_multiple_match.guides, 'Multiple Match Analysis', char(alphabet2(5)));
    xlswrite(batch_filename, {'Strand'}, 'Multiple Match Analysis', char(alphabet1(6)));
    xlswrite(batch_filename, batch_multiple_match.strand, 'Multiple Match Analysis', char(alphabet2(6)));
    xlswrite(batch_filename, {'Chromosome'}, 'Multiple Match Analysis', char(alphabet1(7)));
    xlswrite(batch_filename, batch_multiple_match.chr_id, 'Multiple Match Analysis', char(alphabet2(7)));
    xlswrite(batch_filename, temp_coord_string_start, 'Multiple Match Analysis', char(alphabet1(8)));
    xlswrite(batch_filename, batch_multiple_match.guide_coord_start, 'Multiple Match Analysis', char(alphabet2(8)));
    xlswrite(batch_filename, temp_coord_string_end, 'Multiple Match Analysis', char(alphabet1(9)));
    xlswrite(batch_filename, batch_multiple_match.guide_coord_end, 'Multiple Match Analysis', char(alphabet2(9)));
    xlswrite(batch_filename, temp_coord_string_dsb, 'Multiple Match Analysis', char(alphabet1(10)));
    xlswrite(batch_filename, batch_multiple_match.dsb_coord, 'Multiple Match Analysis', char(alphabet2(10)));
    xlswrite(batch_filename, {'Guide Start (Index)'}, 'Multiple Match Analysis', char(alphabet1(11)));
    xlswrite(batch_filename, batch_multiple_match.guide_index_start, 'Multiple Match Analysis', char(alphabet2(11)));
    xlswrite(batch_filename, {'Guide End (Index)'}, 'Multiple Match Analysis', char(alphabet1(12)));
    xlswrite(batch_filename, batch_multiple_match.guide_index_end, 'Multiple Match Analysis', char(alphabet2(12)));
    xlswrite(batch_filename, {'Double Strand Break (Index)'}, 'Multiple Match Analysis', char(alphabet1(13)));
    xlswrite(batch_filename, batch_multiple_match.dsb_index, 'Multiple Match Analysis', char(alphabet2(13)));
    xlswrite(batch_filename, {'Filename'}, 'Multiple Match Analysis', char(alphabet1(14)));
    xlswrite(batch_filename, batch_multiple_match.filename, 'No Variant Analysis', char(alphabet2(14)));
end

if end_number > starting_number && perform_batch_output == 1 && perform_haplotype_analysis == 1 && output_to_excel == 1
    xlswrite(batch_filename, {'ID'}, 'Haplotype Analysis sgRNA', char(alphabet1(1)));
    xlswrite(batch_filename, batch_unique_haplotype_guides_full.guide_ID_number, 'Haplotype Analysis sgRNA', char(alphabet2(1)));
    xlswrite(batch_filename, {'PAM'}, 'Haplotype Analysis sgRNA', char(alphabet1(2)));
    xlswrite(batch_filename, batch_unique_haplotype_guides_full.pam, 'Haplotype Analysis sgRNA', char(alphabet2(2)));
    xlswrite(batch_filename, {'sgRNA Sequence'}, 'Haplotype Analysis sgRNA', char(alphabet1(3)));
    xlswrite(batch_filename, batch_unique_haplotype_guides_full.guides_only, 'Haplotype Analysis sgRNA', char(alphabet2(3)));
    xlswrite(batch_filename, {'PAM'}, 'Haplotype Analysis sgRNA', char(alphabet1(4)));
    xlswrite(batch_filename, batch_unique_haplotype_guides_full.pam_only, 'Haplotype Analysis sgRNA', char(alphabet2(4)));
    xlswrite(batch_filename, {'sgRNA+PAM'}, 'Haplotype Analysis sgRNA', char(alphabet1(5)));
    xlswrite(batch_filename, batch_unique_haplotype_guides_full.guides, 'Haplotype Analysis sgRNA', char(alphabet2(5)));
    xlswrite(batch_filename, {'Strand'}, 'Haplotype Analysis sgRNA', char(alphabet1(6)));
    xlswrite(batch_filename, batch_unique_haplotype_guides_full.strand, 'Haplotype Analysis sgRNA', char(alphabet2(6)));
    xlswrite(batch_filename, {'Chromosome'}, 'Haplotype Analysis sgRNA', char(alphabet1(7)));
    xlswrite(batch_filename, batch_unique_haplotype_guides_full.chr_id, 'Haplotype Analysis sgRNA', char(alphabet2(7)));
    xlswrite(batch_filename, {'Guide '}, 'Haplotype Analysis sgRNA', char(alphabet1(8)));
    xlswrite(batch_filename, batch_unique_haplotype_guides_full.count, 'Haplotype Analysis sgRNA', char(alphabet2(8)));
    xlswrite(batch_filename, {'Frequency'}, 'Haplotype Analysis sgRNA', char(alphabet1(9)));
    xlswrite(batch_filename, batch_unique_haplotype_guides_full.frequency, 'Haplotype Analysis sgRNA', char(alphabet2(9)));
    xlswrite(batch_filename, {'sgRNA Classification'}, 'Haplotype Analysis sgRNA', char(alphabet1(10)));
    xlswrite(batch_filename, batch_unique_haplotype_guides_full.pam_create, 'Haplotype Analysis sgRNA', char(alphabet2(10)));
    xlswrite(batch_filename, {'N-W Alignment Score'}, 'Haplotype Analysis sgRNA', char(alphabet1(11)));
    xlswrite(batch_filename, batch_unique_haplotype_guides_full.n_mm_wt, 'Haplotype Analysis sgRNA', char(alphabet2(11)));
    xlswrite(batch_filename, {'ID'}, 'Haplotype Analysis sgRNA', char(alphabet1(12)));
    xlswrite(batch_filename, batch_unique_haplotype_guides_full.guide_ID_number_wt, 'Haplotype Analysis sgRNA', char(alphabet2(12)));
    xlswrite(batch_filename, {'Sequence of Associated No Variant sgRNA'}, 'Haplotype Analysis sgRNA', char(alphabet1(13)));
    xlswrite(batch_filename, batch_unique_haplotype_guides_full.guides_only_wt, 'Haplotype Analysis sgRNA', char(alphabet2(13)));
    xlswrite(batch_filename, {'PAM'}, 'Haplotype Analysis sgRNA', char(alphabet1(14)));
    xlswrite(batch_filename, batch_unique_haplotype_guides_full.pam_only_wt, 'Haplotype Analysis sgRNA', char(alphabet2(14)));
    xlswrite(batch_filename, {'sgRNA+PAM'}, 'Haplotype Analysis sgRNA', char(alphabet1(15)));
    xlswrite(batch_filename, batch_unique_haplotype_guides_full.guides_wt, 'Haplotype Analysis sgRNA', char(alphabet2(15)));
    xlswrite(batch_filename, {'Strand'}, 'Haplotype Analysis sgRNA', char(alphabet1(16)));
    xlswrite(batch_filename, batch_unique_haplotype_guides_full.strand_wt, 'Haplotype Analysis sgRNA', char(alphabet2(16)));
    xlswrite(batch_filename, {'Chromosome'}, 'Haplotype Analysis sgRNA', char(alphabet1(17)));
    xlswrite(batch_filename, batch_unique_haplotype_guides_full.chr_id_wt, 'Haplotype Analysis sgRNA', char(alphabet2(17)));
    xlswrite(batch_filename, temp_coord_string_start, 'Haplotype Analysis sgRNA', char(alphabet1(18)));
    xlswrite(batch_filename, batch_unique_haplotype_guides_full.guide_coord_start_wt, 'Haplotype Analysis sgRNA', char(alphabet2(18)));
    xlswrite(batch_filename, temp_coord_string_end, 'Haplotype Analysis sgRNA', char(alphabet1(19)));
    xlswrite(batch_filename, batch_unique_haplotype_guides_full.guide_coord_end_wt, 'Haplotype Analysis sgRNA', char(alphabet2(19)));
    xlswrite(batch_filename, temp_coord_string_dsb, 'Haplotype Analysis sgRNA', char(alphabet1(20)));
    xlswrite(batch_filename, batch_unique_haplotype_guides_full.dsb_coord_wt, 'Haplotype Analysis sgRNA', char(alphabet2(20)));
    xlswrite(batch_filename, {'Guide Start (Index)'}, 'Haplotype Analysis sgRNA', char(alphabet1(21)));
    xlswrite(batch_filename, batch_unique_haplotype_guides_full.guide_index_start_wt, 'Haplotype Analysis sgRNA', char(alphabet2(21)));
    xlswrite(batch_filename, {'Guide End (Index)'}, 'Haplotype Analysis sgRNA', char(alphabet1(22)));
    xlswrite(batch_filename, batch_unique_haplotype_guides_full.guide_index_end_wt, 'Haplotype Analysis sgRNA', char(alphabet2(22)));
    xlswrite(batch_filename, {'Double Strand Break (Index)'}, 'Haplotype Analysis sgRNA', char(alphabet1(23)));
    xlswrite(batch_filename, batch_unique_haplotype_guides_full.dsb_index_wt, 'Haplotype Analysis sgRNA', char(alphabet2(23)));
    xlswrite(batch_filename, {'Filename'}, 'Haplotype Analysis sgRNA', char(alphabet1(24)));
    xlswrite(batch_filename, batch_unique_haplotype_guides_full.filename, 'Haplotype Analysis sgRNA', char(alphabet2(24)));
    
    if end_number > starting_number && perform_batch_output == 1 && perform_haplotype_analysis == 1 && output_to_excel == 1
        num_iterate = max(num_associated_variants_max);
    end
    for h = 1:num_iterate;
        xlswrite(batch_filename, {'RS ID'}, 'Haplotype Analysis sgRNA', char(alphabet1(25+9*(h-1))));
        xlswrite(batch_filename, eval(['batch_unique_haplotype_guides_full',num2str(h),'.rs_id_var']), 'Haplotype Analysis sgRNA', char(alphabet2(25+9*(h-1))));
        xlswrite(batch_filename, {'Ref SNP'}, 'Haplotype Analysis sgRNA', char(alphabet1(26+9*(h-1))));
        xlswrite(batch_filename, eval(['batch_unique_haplotype_guides_full',num2str(h),'.ref_snp_var']), 'Haplotype Analysis sgRNA', char(alphabet2(26+9*(h-1))));
        xlswrite(batch_filename, {'Alt SNP'}, 'Haplotype Analysis sgRNA', char(alphabet1(27+9*(h-1))));
        xlswrite(batch_filename, eval(['batch_unique_haplotype_guides_full',num2str(h),'.alt_snp_var']), 'Haplotype Analysis sgRNA', char(alphabet2(27+9*(h-1))));
        xlswrite(batch_filename, {'Variant Class'}, 'Haplotype Analysis sgRNA', char(alphabet1(28+9*(h-1))));
        xlswrite(batch_filename, eval(['batch_unique_haplotype_guides_full',num2str(h),'.identifier_var']), 'Haplotype Analysis sgRNA', char(alphabet2(28+9*(h-1))));
        xlswrite(batch_filename, {'Variant Chr'}, 'Haplotype Analysis sgRNA', char(alphabet1(29+9*(h-1))));
        xlswrite(batch_filename, eval(['batch_unique_haplotype_guides_full',num2str(h),'.chr_id_var']), 'Haplotype Analysis sgRNA', char(alphabet2(29+9*(h-1))));
        xlswrite(batch_filename, temp_coord_string_start_variant, 'Haplotype Analysis sgRNA', char(alphabet1(30+9*(h-1))));
        xlswrite(batch_filename, eval(['batch_unique_haplotype_guides_full',num2str(h),'.bpstart_var']), 'Haplotype Analysis sgRNA', char(alphabet2(30+9*(h-1))));
        xlswrite(batch_filename, temp_coord_string_end_variant, 'Haplotype Analysis sgRNA', char(alphabet1(31+9*(h-1))));
        xlswrite(batch_filename, eval(['batch_unique_haplotype_guides_full',num2str(h),'.bpend_var']), 'Haplotype Analysis sgRNA', char(alphabet2(31+9*(h-1))));
        xlswrite(batch_filename, {'Variant Start (Index)'}, 'Haplotype Analysis sgRNA', char(alphabet1(32+9*(h-1))));
        xlswrite(batch_filename, eval(['batch_unique_haplotype_guides_full',num2str(h),'.variant_index_start_var']), 'Haplotype Analysis sgRNA', char(alphabet2(32+9*(h-1))));
        xlswrite(batch_filename, {'Variant End (Index)'}, 'Haplotype Analysis sgRNA', char(alphabet1(33+9*(h-1))));
        xlswrite(batch_filename, eval(['batch_unique_haplotype_guides_full',num2str(h),'.variant_index_end_var']), 'Haplotype Analysis sgRNA', char(alphabet2(33+9*(h-1))));
    end
end

if end_number > starting_number && perform_batch_output == 1 && perform_variant_analysis == 1 && output_to_excel == 1 || end_number > starting_number && perform_batch_output == 1 && perform_wgs_analysis == 1 && output_to_excel == 1
    xlswrite(batch_filename, {'ID'}, 'Variant Analysis sgRNA', char(alphabet1(1)));
    xlswrite(batch_filename, batch_unique_variant_guides_full.guide_ID_number, 'Variant Analysis sgRNA', char(alphabet2(1)));
    xlswrite(batch_filename, {'PAM'}, 'Variant Analysis sgRNA', char(alphabet1(2)));
    xlswrite(batch_filename, batch_unique_variant_guides_full.pam, 'Variant Analysis sgRNA', char(alphabet2(2)));
    xlswrite(batch_filename, {'sgRNA Sequence'}, 'Variant Analysis sgRNA', char(alphabet1(3)));
    xlswrite(batch_filename, batch_unique_variant_guides_full.guides_only, 'Variant Analysis sgRNA', char(alphabet2(3)));
    xlswrite(batch_filename, {'PAM'}, 'Variant Analysis sgRNA', char(alphabet1(4)));
    xlswrite(batch_filename, batch_unique_variant_guides_full.pam_only, 'Variant Analysis sgRNA', char(alphabet2(4)));
    xlswrite(batch_filename, {'sgRNA+PAM'}, 'Variant Analysis sgRNA', char(alphabet1(5)));
    xlswrite(batch_filename, batch_unique_variant_guides_full.guides, 'Variant Analysis sgRNA', char(alphabet2(5)));
    xlswrite(batch_filename, {'Strand'}, 'Variant Analysis sgRNA', char(alphabet1(6)));
    xlswrite(batch_filename, batch_unique_variant_guides_full.strand, 'Variant Analysis sgRNA', char(alphabet2(6)));
    xlswrite(batch_filename, {'Chromosome'}, 'Variant Analysis sgRNA', char(alphabet1(7)));
    xlswrite(batch_filename, batch_unique_variant_guides_full.chr_id, 'Variant Analysis sgRNA', char(alphabet2(7)));
    xlswrite(batch_filename, {'Guide '}, 'Variant Analysis sgRNA', char(alphabet1(8)));
    xlswrite(batch_filename, batch_unique_variant_guides_full.count, 'Variant Analysis sgRNA', char(alphabet2(8)));
    xlswrite(batch_filename, {'Frequency'}, 'Variant Analysis sgRNA', char(alphabet1(9)));
    xlswrite(batch_filename, batch_unique_variant_guides_full.frequency, 'Variant Analysis sgRNA', char(alphabet2(9)));
    xlswrite(batch_filename, {'sgRNA Classification'}, 'Variant Analysis sgRNA', char(alphabet1(10)));
    xlswrite(batch_filename, batch_unique_variant_guides_full.pam_create, 'Variant Analysis sgRNA', char(alphabet2(10)));
    xlswrite(batch_filename, {'N-W Alignment Score'}, 'Variant Analysis sgRNA', char(alphabet1(11)));
    xlswrite(batch_filename, batch_unique_variant_guides_full.n_mm_wt, 'Variant Analysis sgRNA', char(alphabet2(11)));
    xlswrite(batch_filename, {'ID'}, 'Variant Analysis sgRNA', char(alphabet1(12)));
    xlswrite(batch_filename, batch_unique_variant_guides_full.guide_ID_number_wt, 'Variant Analysis sgRNA', char(alphabet2(12)));
    xlswrite(batch_filename, {'Sequence of Associated No Variant sgRNA'}, 'Variant Analysis sgRNA', char(alphabet1(13)));
    xlswrite(batch_filename, batch_unique_variant_guides_full.guides_only_wt, 'Variant Analysis sgRNA', char(alphabet2(13)));
    xlswrite(batch_filename, {'PAM'}, 'Variant Analysis sgRNA', char(alphabet1(14)));
    xlswrite(batch_filename, batch_unique_variant_guides_full.pam_only_wt, 'Variant Analysis sgRNA', char(alphabet2(14)));
    xlswrite(batch_filename, {'sgRNA+PAM'}, 'Variant Analysis sgRNA', char(alphabet1(15)));
    xlswrite(batch_filename, batch_unique_variant_guides_full.guides_wt, 'Variant Analysis sgRNA', char(alphabet2(15)));
    xlswrite(batch_filename, {'Strand'}, 'Variant Analysis sgRNA', char(alphabet1(16)));
    xlswrite(batch_filename, batch_unique_variant_guides_full.strand_wt, 'Variant Analysis sgRNA', char(alphabet2(16)));
    xlswrite(batch_filename, {'Chromosome'}, 'Variant Analysis sgRNA', char(alphabet1(17)));
    xlswrite(batch_filename, batch_unique_variant_guides_full.chr_id_wt, 'Variant Analysis sgRNA', char(alphabet2(17)));
    xlswrite(batch_filename, temp_coord_string_start, 'Variant Analysis sgRNA', char(alphabet1(18)));
    xlswrite(batch_filename, batch_unique_variant_guides_full.guide_coord_start_wt, 'Variant Analysis sgRNA', char(alphabet2(18)));
    xlswrite(batch_filename, temp_coord_string_end, 'Variant Analysis sgRNA', char(alphabet1(19)));
    xlswrite(batch_filename, batch_unique_variant_guides_full.guide_coord_end_wt, 'Variant Analysis sgRNA', char(alphabet2(19)));
    xlswrite(batch_filename, temp_coord_string_dsb, 'Variant Analysis sgRNA', char(alphabet1(20)));
    xlswrite(batch_filename, batch_unique_variant_guides_full.dsb_coord_wt, 'Variant Analysis sgRNA', char(alphabet2(20)));
    xlswrite(batch_filename, {'Guide Start (Index)'}, 'Variant Analysis sgRNA', char(alphabet1(21)));
    xlswrite(batch_filename, batch_unique_variant_guides_full.guide_index_start_wt, 'Variant Analysis sgRNA', char(alphabet2(21)));
    xlswrite(batch_filename, {'Guide End (Index)'}, 'Variant Analysis sgRNA', char(alphabet1(22)));
    xlswrite(batch_filename, batch_unique_variant_guides_full.guide_index_end_wt, 'Variant Analysis sgRNA', char(alphabet2(22)));
    xlswrite(batch_filename, {'Double Strand Break (Index)'}, 'Variant Analysis sgRNA', char(alphabet1(23)));
    xlswrite(batch_filename, batch_unique_variant_guides_full.dsb_index_wt, 'Variant Analysis sgRNA', char(alphabet2(23)));
    xlswrite(batch_filename, {'Filename'}, 'Variant Analysis sgRNA', char(alphabet1(24)));
    xlswrite(batch_filename, batch_unique_variant_guides_full.filename, 'Variant Analysis sgRNA', char(alphabet2(24)));
    
    if end_number > starting_number && perform_batch_output == 1 && perform_variant_analysis == 1 && output_to_excel == 1 || end_number > starting_number && perform_batch_output == 1 && perform_wgs_analysis == 1 && output_to_excel == 1
        num_iterate = max(num_associated_variants_max);
    end
    for h = 1:num_iterate;
        xlswrite(batch_filename, {'RS ID'}, 'Variant Analysis sgRNA', char(alphabet1(25+9*(h-1))));
        xlswrite(batch_filename, eval(['batch_unique_variant_guides_full',num2str(h),'.rs_id_var']), 'Variant Analysis sgRNA', char(alphabet2(25+9*(h-1))));
        xlswrite(batch_filename, {'Ref SNP'}, 'Variant Analysis sgRNA', char(alphabet1(26+9*(h-1))));
        xlswrite(batch_filename, eval(['batch_unique_variant_guides_full',num2str(h),'.ref_snp_var']), 'Variant Analysis sgRNA', char(alphabet2(26+9*(h-1))));
        xlswrite(batch_filename, {'Alt SNP'}, 'Variant Analysis sgRNA', char(alphabet1(27+9*(h-1))));
        xlswrite(batch_filename, eval(['batch_unique_variant_guides_full',num2str(h),'.alt_snp_var']), 'Variant Analysis sgRNA', char(alphabet2(27+9*(h-1))));
        xlswrite(batch_filename, {'Variant Class'}, 'Variant Analysis sgRNA', char(alphabet1(28+9*(h-1))));
        xlswrite(batch_filename, eval(['batch_unique_variant_guides_full',num2str(h),'.identifier_var']), 'Variant Analysis sgRNA', char(alphabet2(28+9*(h-1))));
        xlswrite(batch_filename, {'Variant Chr'}, 'Variant Analysis sgRNA', char(alphabet1(29+9*(h-1))));
        xlswrite(batch_filename, eval(['batch_unique_variant_guides_full',num2str(h),'.chr_id_var']), 'Variant Analysis sgRNA', char(alphabet2(29+9*(h-1))));
        xlswrite(batch_filename, temp_coord_string_start_variant, 'Variant Analysis sgRNA', char(alphabet1(30+9*(h-1))));
        xlswrite(batch_filename, eval(['batch_unique_variant_guides_full',num2str(h),'.bpstart_var']), 'Variant Analysis sgRNA', char(alphabet2(30+9*(h-1))));
        xlswrite(batch_filename, temp_coord_string_end_variant, 'Variant Analysis sgRNA', char(alphabet1(31+9*(h-1))));
        xlswrite(batch_filename, eval(['batch_unique_variant_guides_full',num2str(h),'.bpend_var']), 'Variant Analysis sgRNA', char(alphabet2(31+9*(h-1))));
        xlswrite(batch_filename, {'Variant Start (Index)'}, 'Variant Analysis sgRNA', char(alphabet1(32+9*(h-1))));
        xlswrite(batch_filename, eval(['batch_unique_variant_guides_full',num2str(h),'.variant_index_start_var']), 'Variant Analysis sgRNA', char(alphabet2(32+9*(h-1))));
        xlswrite(batch_filename, {'Variant End (Index)'}, 'Variant Analysis sgRNA', char(alphabet1(33+9*(h-1))));
        xlswrite(batch_filename, eval(['batch_unique_variant_guides_full',num2str(h),'.variant_index_end_var']), 'Variant Analysis sgRNA', char(alphabet2(33+9*(h-1))));
    end
end

disp('DNA Striker analysis complete!')

