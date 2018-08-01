#!/bin/env julia

using CLI, JLD, Distributions, HypothesisTests

immutable Region
	start::Int32
	stop::Int32
end

# Print comment lines to stdout and return
# first non-comment line
function skip_to_header(vcf_file)
	while true
		line = readline(vcf_file); 
		if !startswith( line, "#") return line; end		
		println(line);
	end
end

drop_chr_prefix(chr) = startswith(chr, "chr") ? chr[4:end] : chr

is_protein_altering(effect) = ismatch(
	r"Missense|[Ff]rameshift|Stopgain|Stoploss|Splice", effect);

function mutation_class(ref, alt)
	if ref == "C"
		if alt == "A"; return 1; end
		if alt == "G"; return 2; end
		if alt == "T"; return 3; end
	elseif ref == "G"
		if alt == "T"; return 1; end
		if alt == "C"; return 2; end
		if alt == "A"; return 3; end
	elseif ref == "T"
		if alt == "A"; return 4; end
		if alt == "C"; return 5; end
		if alt == "G"; return 6; end
	elseif ref == "A"
		if alt == "T"; return 4; end
		if alt == "G"; return 5; end
		if alt == "C"; return 6; end
	end

	if length(ref) > 1 || ref == "-" || length(alt) > 1 || alt == "-"
		return 0
	else
		error("Mutation class could not be defined.")
	end
end




function somatic(vcf_file::IO, samples_file::IO; alt_reads=3, alt_frac=0.1, 
	test_ref_ratio=5, test_bg_ratio=20, test_bg_alpha=1.0, ref_reads=0,
	keep_all=false, min_mapq=10, min_sidedness=6, keep_old_stars::Bool=false)

	header = skip_to_header(vcf_file); println(header)
	header_cols = split(header, '\t') 
	notes_col = findfirst(header_cols .== "NOTES")
	alt_col = findfirst(header_cols .== "ALT")
	ref_col = findfirst(header_cols .== "REF")

	samples = header_cols[notes_col+1:end]; 
	n_samples = length(samples);

	# Read samples file and separate the row header
	d = readdlm(samples_file, '\t'); 
	sample_file_header = d[1, :][:]; 
	d = d[2:end, :];

	sample_test_col = findfirst(sample_file_header .== "TEST")
	sample_ref_col = findfirst(sample_file_header .== "REF")
	assert(sample_test_col != 0 && sample_ref_col != 0)

	# Check that all samples named in the design file are present in the VCF
	for s in unique(d[:, [sample_test_col, sample_ref_col]])
		if s != "" && !(s in samples)
			error("Sample $s is listed in the design file but missing from VCF.")
		end
	end

	# Background error is estimated based on paired *and* unpaired reference
	# samples.
	ref_samples = filter(s -> s != "", unique(d[:, sample_ref_col]))
	background = indexin(ref_samples, samples)

	# Get indices of test and ref sample rows in the samples file
	test_rows = d[:, sample_test_col] .!= ""
	test = indexin(d[test_rows, sample_test_col], samples)
	n_tests = length(test);
	ref = indexin(d[test_rows, sample_ref_col], samples)
	info("Analyzing $n_tests tumor-normal pairs...")

	if !isempty(findin(background, test))
		warn("Some background reference samples also used as test samples.")
	end
	info("Calculating background error rate based on $(length(background)) samples...")

	# Read alt_frac from design file or command line parameter.
	alt_frac_col = findfirst(sample_file_header .== "ALT_FRAC")
	if alt_frac_col == 0
		info("Using a mutant allele fraction threshold of $(alt_frac * 100)% for all samples...")
		alt_frac = ones(n_tests) * alt_frac
	else
		alt_frac = convert(Array{Float64}, d[test_rows, alt_frac_col])
	end

	# TODO: Allow sample-specific values to be specified in design file for
	# alt_reads, test_ref_ratio and test_bg_ratio as well.

	for line in eachline(vcf_file)

		cols = split(line, '\t')
		alt = zeros(Int32, n_samples); 
		total = zeros(Int32, n_samples);
		# CALL2 support
		mapq = zeros(Int32, n_samples);
		baseq = zeros(Int32, n_samples);
		sidedness = zeros(Int32, n_samples);

		call2_format = false
		is_indel = length(cols[alt_col]) != length(cols[ref_col])

		for s in 1:n_samples
			col = cols[notes_col + s]			
			gt = split(cols[notes_col + s], ':')
			alt[s] = parse(Int32, gt[1])
			total[s] = parse(Int32, gt[2])
			# CALL2 support			
			if length(gt) >= 5
				mapq[s] = parse(Int32, gt[3])
				baseq[s] = parse(Int32, gt[4])
				sidedness[s] = parse(Int32, gt[5])
				call2_format = true
			end
		end

		fracs = alt ./ total
		bg_alt = sum(alt[background])
		bg_total = sum(total[background])
		avg_ref_frac = bg_alt / bg_total

		somatic = falses(n_samples)
		for k in 1:n_tests
			t = test[k]; # Index of test sample col
			r = ref[k];  # Index of ref sample col

			# Perform checks to determine somatic status of mutation
			if fracs[t] < alt_frac[k]; continue; end
			if r != 0 && fracs[t] < test_ref_ratio * fracs[r]; continue; end
			if fracs[t] < test_bg_ratio * avg_ref_frac; continue; end
			if alt[t] < alt_reads; continue; end
			if r != 0 && total[r] < ref_reads; continue; end
			if call2_format
				if !is_indel && mapq[t] < min_mapq; continue; end 
				if sidedness[t] < min_sidedness; continue; end 
			end

			#if test_bg_alpha < 1.0
			#	if pvalue(ChisqTest([alt[t] total[t] - alt[t]; bg_alt bg_total - bg_alt])) > test_ref_alpha; continue; end
			#end
			#chisq_p = pvalue(ChisqTest(
			#	[alt[t] total[t] - alt[t]; bg_alt bg_total - bg_alt]))
			#baseq[t] = round(Int32, min(-10 * log10(chisq_p), 100))

			somatic[t] = true
		end

		if keep_old_stars
			somatic |= map(c -> endswith(c, '*'), cols[notes_col+1:end])
		end

		if !any(somatic) && keep_all == false; continue; end
		
		print(join(cols[1:notes_col], '\t'))
		@printf("BG: %.3f%%. ", avg_ref_frac * 100)
		for s in 1:n_samples
			if call2_format
				print("\t$(alt[s]):$(total[s]):$(mapq[s]):$(baseq[s]):$(sidedness[s])")				
			else 
				print("\t$(alt[s]):$(total[s])")
			end

			if somatic[s]; print(":*"); end
		end
		println()
	end
end

function germline(vcf_file::IO, ref_regexps...; alt_reads=5, alt_frac=0.15, bg_ratio=20)
	line = skip_to_header(vcf_file); println(line)
	sample_col = findfirst(split(line, '\t') .== "NOTES") + 1
	samples = split(line, '\t')[sample_col:end]
	germline = falses(samples)
	for regex in ref_regexps; germline .|= ismatch.(Regex(regex), samples); end
	for line in eachline(vcf_file)
		cols = split(line, '\t')
		gtypes = [split(c, ':') for c in cols[sample_col:end]]
		total = [parse(Int32, gt[2]) for gt in gtypes]
		alt = [parse(Int32, gt[1]) for gt in gtypes]
		fracs = alt ./ total
		alt_gt = (fracs .>= alt_frac) .& (alt .>= alt_reads) .& germline
		if !any(alt_gt); continue; end

		# There might be non-reference samples that are not in alt_gt because
		# of the alt_reads threshold. We don't want to include these
		# samples in the background, so we only use alt_frac when picking
		# background samples.
		avg_alt = sum(alt[alt_gt]) / sum(total[alt_gt])
		bg = (fracs .< alt_frac) .& germline
		avg_ref = sum(alt[bg]) / sum(total[bg])
		if avg_alt / avg_ref < bg_ratio; continue; end

		print(join(cols[1:sample_col-1], '\t'))
		for s in 1:length(total)
			print("\t$(alt[s]):$(total[s])")
			if alt_gt[s]; print(":*"); end
		end
		print('\n')
	end
end

function above_background(vcf_file::IO, ref_regex...; alt_reads=3, alt_frac=0.1, bg_ratio=20)
	line = skip_to_header(vcf_file); println(line)
	sample_col = findfirst(split(line, '\t') .== "NOTES") + 1
	samples = split(line, '\t')[sample_col:end]

	bg = falses(samples)
	for regex in ref_regex; bg .|= ismatch.(Regex(regex), samples); end
	info("Calculating background based on $(sum(bg)) samples...")

	for line in eachline(vcf_file)
		cols = split(line, '\t')
		gtypes = map(c -> split(c, ':'), cols[sample_col:end])
		total = map(gt -> parse(Int, gt[2]), gtypes)
		alt = map(gt -> parse(Int, gt[1]), gtypes)
		fracs = alt ./ total
		avg_bg_frac = sum(alt[bg]) / sum(total[bg])
		alt_gt = (alt .>= alt_reads) .&
			(fracs .>= max(alt_frac, avg_bg_frac * bg_ratio))
		if !any(alt_gt); continue; end

		print(join(cols[1:sample_col-1], '\t'))
		for s in 1:length(total)
			print("\t$(alt[s]):$(total[s])")
			if alt_gt[s]; print(":*"); end
		end
		println()
	end
end

function heterozygous_snps(vcf_file::IO, ref_regex...; min_depth=30, min_frac=0.3, bg_fraction=0.1)
	line = skip_to_header(vcf_file); println(line)
	sample_col = findfirst(split(line, '\t') .== "NOTES") + 1
	samples = split(line, '\t')[sample_col:end]; S = length(samples)

	ref = falses(samples)
	for regex in ref_regex; ref .|= ismatch.(Regex(regex), samples); end
	info("Searching for heterozygous SNPs in $(sum(ref)) germline samples...")

	for line in eachline(vcf_file)
		cols = split(line, '\t')
		alt = zeros(Int, S); total = zeros(Int, S); starred = falses(S);
		unclear = 0
		for (s, col) in enumerate(cols[sample_col:end])
			parts = split(col, ':')
			alt[s] = parse(Int, parts[1])
			total[s] = parse(Int, parts[2])
			if ref[s] == false; continue; end
			if total[s] < min_depth; continue; end
			
			if alt[s] / total[s] < min_frac
				error_lower = Binomial(total[s], 0.01)
				if 1 - cdf(error_lower, alt[s] - 1) <= 0.05; unclear += 1; end
			elseif alt[s] / total[s] > 1 - min_frac
				error_upper = Binomial(total[s], 0.99)
				if cdf(error_upper, alt[s]) <= 0.05; unclear += 1; end
			end

			if unclear / S >= bg_fraction
				starred[:] = false; break;
			end

			hetz = Binomial(total[s], 0.5)
			if 0.05 <= cdf(hetz, alt[s]) <= 0.95; starred[s] = true; end
		end

		if !any(starred); continue; end

		print(join(cols[1:sample_col-1], '\t'))
		for s in 1:S
			print("\t$(alt[s]):$(total[s])")
			if starred[s]; print(":*"); end
		end
		print('\n')
	end
end

annovar_valid_funcs = Set(["exonic", "splicing", "intergenic", "intronic",
	"upstream", "downstream", "UTR5", "UTR3", "ncRNA_exonic", "ncRNA_intronic",
	"ncRNA_splicing", ""])

function translate_annovar_effect(gene, func, exonic_func, aa_change)
	for f in split(func, ';')
		if !in(f, annovar_valid_funcs)
			warn("Unrecognized variant effect '$f'.")
		end
	end

	effects = []
	if contains(func, "splicing"); push!(effects, "Splice site"); end
	if contains(func, "exonic")
		details = unique(map(m -> m.captures[1],
			eachmatch(r":(p\..+?)(,|$)", aa_change)))
		if exonic_func == "synonymous SNV"
			push!(effects, "Synonymous $(join(details, ", "))")
		elseif exonic_func == "nonsynonymous SNV"
			push!(effects, "Missense $(join(details, ", "))")
		elseif exonic_func == "stopgain"
			push!(effects, "Stopgain $(join(details, ", "))")
		elseif exonic_func == "stoploss"
			push!(effects, "Stoploss $(join(details, ", "))")
		elseif contains(exonic_func, "nonframeshift")
			push!(effects, "Non-frameshift indel $(join(details, ", "))")
		elseif contains(exonic_func, "frameshift")
			push!(effects, "Frameshift $(join(details, ", "))")
		elseif contains(func, "ncRNA_exonic")
			push!(effects, "Exonic (ncRNA)")
		elseif exonic_func == "unknown"
			push!(effects, "Exonic (unknown)")
		else
			error("Unrecognized effect: $([func, exonic_func])")
		end
	end
	if contains(func, "UTR3"); push!(effects, "3'-UTR"); end
	if contains(func, "UTR5"); push!(effects, "5'-UTR"); end
	if contains(func, "upstream"); push!(effects, "Upstream"); end
	if contains(func, "downstream"); push!(effects, "Downstream"); end
	if contains(func, "intronic"); push!(effects, "Intronic"); end
	if func == "intergenic"; push!(effects, "Intergenic"); gene = ""; end

	return (gene, join(effects, ". "))
end

function predict_effect(vcf_file::IO; genome="~/tools/annovar-2016-02-01/humandb/hg38")
	genome_version = basename(genome)
	humandb_dir = expanduser(dirname(genome))
	tmp_prefix = "annovar-$(hex(rand(0:2^16-1), 4))"
	original_path = "$(tmp_prefix).original"
	reformatted_path = "$(tmp_prefix).avinput"

	# Reformat the input file to an ANNOVAR-compatible format. In particular,
	# indels must be changed so they use '-' as the ref or alt allele.
	# Since the input file could be a pipe, we must write the original VCF
	# into a separate file that will be used to build the final annotated VCF.
	original_file = open(original_path, "w")
	reformatted_file = open(reformatted_path, "w")	
	headers = split(skip_to_header(vcf_file), '\t')
	
	for line in eachline(vcf_file)
		
		write(original_file, line, '\n')
		c = split(line, '\t')
		start = parse(Int, c[2]); stop = start; ref = c[3]; alt = c[4]

		if length(ref) == 1 && length(alt) > 1 && ref[1] == alt[1]
			# Simple insertion
			ref = "-"; alt = alt[2:end];
		elseif length(ref) > 1 && length(alt) == 1 && ref[1] == alt[1]
			# Simple deletion
			ref = ref[2:end]; alt = "-"; start += 1;
		elseif length(ref) > 1 && length(alt) > 1
			# Block substitution
		end
		stop = start + length(ref) - 1
		write(reformatted_file, "$(c[1])\t$(start)\t$(stop)\t$(ref)\t$(alt)\n")
	end
	close(original_file); close(reformatted_file); close(vcf_file)

	run(pipeline(`table_annovar.pl $(reformatted_path) $(humandb_dir) --buildver $(genome_version) --remove --outfile $(tmp_prefix) --operation g --protocol refGene`, stderr=DevNull))

	anno_path = "$(tmp_prefix).$(genome_version)_multianno.txt"
	original_file = open(original_path);
	anno_file = open(anno_path); readline(anno_file)
	print("CHROM\tPOSITION\tREF\tALT\tGENE\tEFFECT\tNOTES\t")
	println(join(headers[6:end], '\t'))
	for line in eachline(anno_file)
		c = split(line, '\t')
		if length(c) >= 7
			gene, effect = translate_annovar_effect(c[7], c[6],
				length(c) >= 9 ? c[9] : "", length(c) >= 10 ? c[10] : "")
		else
			gene = ""; effect = "";
		end

		orig_line = readline(original_file)
		tab = search(orig_line, '\t')
		for k in 1:3; tab = search(orig_line, '\t', tab + 1); end
		print(orig_line[1:tab])
		print("$(gene)\t$(effect)\t")
		println(orig_line[tab+1:end])
	end
	close(anno_file); rm(original_path); rm(reformatted_path); rm(anno_path);
end

type TextAnnotations
	position::Vector{Int32}
	change::Vector{Int16}
	annotation::Vector{String}
	TextAnnotations() = new([], [], [])
end

type NumericAnnotations
	position::Vector{Int32}
	change::Vector{Int16}
	annotation::Vector{Float32}
	NumericAnnotations() = new([], [], [])
end

base_num = Dict{String, Int16}("A" => 1, "C" => 2, "G" => 3, "T" => 4)
encode_ref_alt(ref, alt) = get(base_num, ref, 0) * 16 + get(base_num, alt, 0)

function build_annotation_database(anno_path, db_path; prefix="",
	format="text")
	db = jldopen(db_path, "w", compress=true)
	anno_file = zopen(anno_path)
	if format == "text"
		build_annotation_database_text(anno_file, db)
	elseif format == "numeric"
		build_annotation_database_numeric(anno_file, db)
	elseif format == "blank"
		build_annotation_database_blank(anno_file, db)
	else
		error("Unrecognized annotation database format '$format' requested.")
	end
	write(db, "prefix", prefix)
	write(db, "format", format)
	close(db)
end

function build_annotation_database_blank(anno_file, db)
	chr = ""; anno = TextAnnotations()
	for line in eachline(anno_file)
		c = split(rstrip(line, '\t'))
		if c[1] != chr
			if chr != ""
				assert(issorted(anno.position))
				write(db, chr, anno)
			end
			anno = TextAnnotations(); chr = String(c[1])
			println(STDERR, "Chromosome $(chr)...")
		end
		push!(anno.position, parse(Int32, c[2]))
		push!(anno.change, encode_ref_alt(c[3], c[4]))
	end
	write(db, chr, anno)
end

function build_annotation_database_text(anno_file, db)
	chr = ""; anno = TextAnnotations()
	for line in eachline(anno_file)
		c = split(rstrip(line, '\t'))
		if c[1] != chr
			if chr != ""
				assert(issorted(anno.position))
				write(db, chr, anno)
			end
			anno = TextAnnotations(); chr = String(c[1])
			println(STDERR, "Chromosome $(chr)...")
		end
		push!(anno.position, parse(Int32, c[2]))
		push!(anno.change, encode_ref_alt(c[3], c[4]))
		push!(anno.annotation, c[5])
	end
	write(db, chr, anno)
end

function build_annotation_database_numeric(anno_file, db)
	chr = ""; anno = NumericAnnotations()
	for line in eachline(anno_file)
		c = split(rstrip(line, '\t'))
		if c[1] != chr
			if chr != ""
				assert(issorted(anno.position))
				write(db, chr, anno)
			end
			anno = NumericAnnotations(); chr = String(c[1])
			println(STDERR, "Chromosome $(chr)...")
		end
		push!(anno.position, parse(Int32, c[2]))
		push!(anno.change, encode_ref_alt(c[3], c[4]))
		push!(anno.annotation, parse(Float32, c[5]))
	end
	write(db, chr, anno)
end

function row_to_key( row::String)	
	# Firstly, compare non-numeric parts (case-insensitive)
	# Next, compare numeric parts as integer
	# Next, compare position as integer
	chrom, pos = split(row, '\t')[1:2]
	return (filter(x -> !isnumber(x), lowercase(chrom)),
		    parse(Int32, filter(x -> isnumber(x), '0'*chrom)),
		    parse(Int32, pos))
end

# Sorts and prints VCF file rows by their chromosomal position
# Does not allow comment lines anywhere else than the beginning
function sort_by_position(vcf_file::IO)

	println( skip_to_header(vcf_file));
	rows = readlines(vcf_file) 

	sort!(rows, by=x -> row_to_key(x))
	for row in rows println(row); end
end



function annotate(vcf_file::IO, db_path)
	db = jldopen(db_path, "r")
	prefix = read(db, "prefix")
	format = read(db, "format")
	headers = split(skip_to_header(vcf_file), '\t')
	notes_col = findfirst(headers .== "NOTES")
	println(join(headers, '\t'))

	chr = ""; anno = nothing
	for line in eachline(vcf_file)

		c = split(line, '\t')
		if drop_chr_prefix(c[1]) != chr
			chr = String(drop_chr_prefix(c[1]))
			if any(x -> chr == x, names(db))   # in() does not work
				anno = read(db, chr)
			else
				warn("No annotations found for chromosome $(chr).")
				anno = NumericAnnotations()
			end
		end
		matches = searchsorted(anno.position, parse(Int32, c[2]))
		if !isempty(matches)
			change = encode_ref_alt(c[3], c[4])
			for m in matches
				if anno.change[m] == change
					if format == "numeric"
						# See https://github.com/JuliaLang/julia/issues/14331
						c[notes_col] *= rstrip(@sprintf("%s:%g", prefix,
							anno.annotation[m])) * ". "
					elseif format == "text"
						c[notes_col] *= "$(prefix):$(anno.annotation[m]). "
					else
						c[notes_col] *= "$(prefix). "
					end
				end
			end
		end
		println(join(c, '\t'))
	end
end

# Filter VCF file lines based on if they are protein altering or not
# ARG: invert: Set to false to output protein altering mutations"
#              Set to true to output non protein altering mutations"
function protein_altering(vcf_file::IO; invert=false)
	line = skip_to_header(vcf_file); println(line)
	func_col = findfirst(split(line, '\t') .== "EFFECT")
	for line in eachline(vcf_file)
		cols = split(line, '\t')
		alters = ismatch(r"Missense|rameshift|Stopgain|Stoploss|Splice",
			cols[func_col])
		if !alters == invert; println(line); end
	end
end

function select_samples(invert::Bool, vcf_file::IO, regexps...)
	headers = split(skip_to_header(vcf_file), '\t')
	sample_col = findfirst(headers .== "NOTES") + 1
	keep = falses(length(headers))
	keep[1:sample_col-1] = true
	for col in sample_col:length(headers)
		keep[col] = any(rx -> ismatch(Regex(rx), headers[col]), regexps)
	end
	if invert; keep[sample_col:end] = .!keep[sample_col:end]; end
	println(join(headers[keep], '\t'))
	for line in eachline(vcf_file)
		cols = split(line, '\t')
		println(join(cols[keep], '\t'))
	end
end
keep_samples(vcf_file::IO, regexps...) = select_samples(false, vcf_file, regexps...)
discard_samples(vcf_file::IO, regexps...) = select_samples(true, vcf_file, regexps...)

function discard_shallow(vcf_file::IO, min_median_depth::Float64)
	line = skip_to_header(vcf_file); println(line);
	sample_col = findfirst(split(line, '\t') .== "NOTES") + 1
	for line in eachline(vcf_file)
		total = map(c -> parse(Int, split(c, ':')[2]),
			split(line, '\t')[sample_col:end])
		if median(total) >= min_median_depth; println(line); end
	end
end

function discard_blacklisted(vcf_file::IO, blacklist::IO)
	d = readdlm(blacklist, '\t', String)
	blacklist = Set(map(r -> d[r, 1:4][:], 1:size(d, 1)))

	line = skip_to_header(vcf_file); println(line);
	for line in eachline(vcf_file)
		c = split(line, '\t')
		if in(c[1:4], blacklist); continue; end
		println(line)
	end
end

function discard_indels(vcf_file::IO)
	line = skip_to_header(vcf_file); println(line);
	for line in eachline(vcf_file)
		c = split(line, '\t')
		if length(c[3]) == 1 && length(c[4]) == 1 && c[3] != "-" && c[4] != "-"
			println(line)
		end
	end
end

function discard_if_frequency_above(vcf_file::IO, db_path, threshold::Float64)
	db = jldopen(db_path, "r")
	println(skip_to_header(vcf_file))
	chr = ""; anno = nothing
	for line in eachline(vcf_file)
		c = split(line, '\t')
		if drop_chr_prefix(c[1]) != chr
			chr = String(drop_chr_prefix(c[1])); anno = read(db, chr)
		end
		matches = searchsorted(anno.position, parse(Int32, c[2]))
		if !isempty(matches)
			change = encode_ref_alt(c[3], c[4])
			if any(m -> anno.change[m] == change &&
				anno.annotation[m] >= threshold, matches)
				continue;
			end
		end
		println(join(c, '\t'))
	end
end

function inside(vcf_file::IO, region)
	m = match(r"([a-zA-Z0-9]+):(\d+)-(\d+)", region)
	if m == nothing
		if isfile(region) || isfifo(region)
			inside_bed_regions(vcf_file, region)
			return
		else
			error("Invalid region specified.")
		end
	end

	chr = m.captures[1]
	start = parse(Int, m.captures[2])
	stop = parse(Int, m.captures[3])

	line = skip_to_header(vcf_file); println(line)
	sample_col = findfirst(split(line, '\t') .== "NOTES") + 1
	for line in eachline(vcf_file)
		c = split(line, '\t'); pos = parse(Int, c[2])
		if c[1] == chr && start <= pos <= stop; println(line); end
	end
end

function inside_bed_regions(vcf_file::IO, regions_path)
	regions = Dict{String, Array{Region}}()
	for line in eachline(open(regions_path))
		c = split(rstrip(line), '\t')
		push!(get!(() -> Array{Region}(0), regions, c[1]),
			Region(parse(Int, c[2])+1, parse(Int, c[3])))
	end

	line = skip_to_header(vcf_file); println(line)
	sample_col = findfirst(split(line, '\t') .== "NOTES") + 1 + 1
	for line in eachline(vcf_file)
		c = split(line, '\t'); pos = parse(Int, c[2])
		chr_regions = get(() -> Array{Region}(0), regions, c[1])
		if any(r -> r.start <= pos <= r.stop, chr_regions)
			println(line)
		end
	end
end

function discard_contingent(vcf_path, alpha::Float64)
	vcf_file = zopen(vcf_path)
	line = skip_to_header(vcf_file); print(line)
	sample_col = findfirst(split(line, '\t') .== "NOTES") + 1
	for line in eachline(vcf_file)
		gtypes = map(c -> split(c, ':'), split(line, '\t')[sample_col:end])
		total = map(gt -> parse(Int, gt[2]), gtypes)
		alt = map(gt -> parse(Int, gt[1]), gtypes)
		valid = total .>= 1
		p = any(valid) ? pvalue(ChisqTest([alt[valid] total[valid] - alt[valid]])) : 1
		if alpha == 1
			print(line[1:end-1]); print("\t$(p)\n")
		elseif p <= alpha
			print(line)
		end
	end
end

# Adds a warning to the Notes column if the mappability in the reference
# genome is below a specified treshold (-> low mappability)
function mappability(vcf_file::IO, track_path::String, threshold::Float64)
	tracks = jldopen(track_path, "r");
	headers = split(skip_to_header(vcf_file), '\t')
	notes_col = findfirst(headers .== "NOTES")
	println(join(headers, '\t'))

	# We assume that the VCF file is sorted by chromosome.
	curr_chr = ""; track = zeros(UInt8, 0);
	for line in eachline(vcf_file)
		c = split(line, '\t')
		if c[1] != curr_chr
			curr_chr = String(c[1])
			track = read(tracks, curr_chr)
		end
		mappability = track[parse(Int, c[2])]
		if mappability < threshold
			c[notes_col] = "$(c[notes_col])Mappability < $(threshold)%. "
		end
		println(join(c, '\t'))
	end
end

# Adds a warning to the NOTES column if there are indels within
# dist distance of the alt read (line in VCF file).
#
# ARG: dist : largest distance that is considered
#             being "nearby"
function nearby_indels(vcf_file::IO; dist=10)

	line = skip_to_header(vcf_file);
	println(line)
	headers = split(line, '\t')
	notes_col = findfirst(split(line, '\t') .== "NOTES")
	warn_msg = "INDEL within $dist bases. "

	#notes_col = findfirst(split(skip_to_header(vcf_file), '\t') .== "NOTES")

	prev_lines = Array{Array{String}}(0)

	for line in eachline(vcf_file)

		if startswith(line, "#") println(line); continue; end 
		c = split(line, '\t')
		pos = parse(Int, c[2])
		starred = map(x -> endswith(x, '*'), c[notes_col+1:end])

		# Forget old lines that are so far behind they cannot be affected.
		keep = trues(length(prev_lines))
		for (k, prev_line) in enumerate(prev_lines)
			# On different chromosome or too distant
			if c[1] != prev_line[1] || parse(Int, prev_line[2]) < pos - dist
				keep[k] = false
				println(join(prev_line, '\t'))
			end
		end
		prev_lines = prev_lines[keep]

		# Add warning to old lines if current line is an indel.
		if length(c[3]) != length(c[4])
			for prev_line in prev_lines
				if any(s -> starred[s] && endswith(prev_line[notes_col+s], '*'), 1:length(starred))
					prev_line[notes_col] *= warn_msg
				end
			end
		end

		# Add warning to current line if any old lines were indels.
		for prev_line in prev_lines
			if length(prev_line[3]) != length(prev_line[4]) &&
				any(s -> starred[s] && endswith(prev_line[notes_col+s], '*'), 1:length(starred))
				c[notes_col] *= warn_msg
				break
			end
		end

		push!(prev_lines, c)
	end

	for prev_line in prev_lines; println(join(prev_line, '\t')); end
end

function statistics(vcf_path)
	vcf_file = zopen(vcf_path)
	headers = split(rstrip(skip_to_header(vcf_file)), '\t')
	sample_col = findfirst(headers .== "NOTES") + 1
	samples = headers[sample_col:end]; S = length(samples)
	total_subs = zeros(Int, S)
	possible_subs = zeros(Int, S, 6)
	total_indels = zeros(Int, S)

	for line in eachline(vcf_file)
		c = split(rstrip(line), '\t')
		alt = map(gt -> endswith(gt, '*'), c[sample_col:end])
		mut_class = mutation_class(c[3], c[4])
		if mut_class > 0
			total_subs += alt
			possible_subs[:, mut_class] += alt
		elseif mut_class == 0
			total_indels += alt
		end
	end

	println("SAMPLE\tTOTAL SUBSTITUTIONS\tTOTAL INDELS\tC>A\tC>G\tC>T\tT>A\tT>C\tT>G")
	for s in 1:S
		print("$(samples[s])\t$(total_subs[s])\t$(total_indels[s])")
		for k in 1:6; print("\t$(possible_subs[s, k])"); end
		println()
	end
end


# CALL2 output format: "alt_reads:total_reads:mapq:baseq:sidedness:significant"
function fractions(vcf_file::IO)

	line = skip_to_header(vcf_file); println(line)
	sample_col = findfirst(split(line, '\t') .== "NOTES") + 1
	for line in eachline(vcf_file)
		cols = split(line, '\t')
		gtypes = [split(c, ':') for c in cols][sample_col:end]

		total = [parse(Int32, gt[2]) for gt in gtypes]
		alt = [parse(Int32, gt[1]) for gt in gtypes]

		print(join(cols[1:(sample_col-1)], '\t'))

		for s in 1:length(gtypes)
			frac = alt[s] / total[s]			
			@printf("\t%.1f%% (%d)", isnan(frac) ? 0.0 : frac * 100, total[s])
			#CALL2 support
			if length(gtypes[s]) >= 5 && parse( Int32, gtypes[s][ 3]) > 0; @printf(" [mq:%s,bq:%s,sd:%s]", gtypes[s][3], gtypes[s][4], gtypes[s][5] ); end						
			if length(gtypes[s]) >= 3 && gtypes[s][ length(gtypes[s])] == "*"; print(" *"); end
		end
		println()
	end
end

function error_rate(bam_path, genome_path; min_mapq=10)
	bam_path = expanduser(bam_path)
	genome_path = expanduser(genome_path)
	spileup = expanduser("~/tools/pypette/compiled/spileup")

	# 2D histogram, mutant read count on rows, mutant read fraction on columns
	# Histogram bins are left-exclusive, right-inclusive.
	read_thresholds = 0:100
	frac_thresholds = 0:0.001:0.10
	histogram = zeros(Int, length(read_thresholds), length(frac_thresholds))

	for line in eachline(pipeline(`samtools mpileup -d 100000 -A -x -R -sB -q$(min_mapq) -f $(genome_path) $(bam_path)`, `$(spileup) 0 $(min_mapq)`))
		cols = split(rstrip(line), '\t')
		if length(cols) < 4; continue; end
		#if cols[3] == "N"; print(line); continue; end
		assert(length(cols) == 4)
		assert(cols[3] != "N")

		pileup = cols[4]
		if pileup == ""; continue; end
		t = split(pileup, ' ')

		total_reads = 0
		alt_reads = 0

		for a in 1:3:length(t)
		 	count = parse(Int, t[a+1])
		 	total_reads += count
		 	if t[a] != "."; alt_reads += count; end
		end

		alt_frac = alt_reads / total_reads
		read_bin = Int(min(alt_reads + 1, 101))
		frac_bin = min(ceil(Int, alt_frac / 0.001) + 1, 101)
		histogram[read_bin, frac_bin] += 1

		# allele_reads = zeros(Int32, max_alleles)
		# alleles = Array{String}(0)

		# pileup = cols[4]
		# if pileup == ""; continue; end
		# t = split(pileup, ' ')
		# for a in 1:3:length(t)
		# 	count = int(t[a+1])
		# 	total_reads[s] += count
		# 	if t[a] == "."; continue; end
		# 	idx = findfirst(allele -> allele == t[a], alleles)
		# 	if idx == 0 && length(alleles) >= max_alleles; continue; end
		# 	if idx == 0; push!(alleles, t[a]); idx = length(alleles); end
		# 	allele_reads[s, idx] = count
		# end

		# for a in 1:length(alleles)
		# end
	end

	writedlm(STDOUT, histogram)
end

function group_by_snp_profile(hetz_snp_vcf::IO)
	headers = split(rstrip(skip_to_header(hetz_snp_vcf)), '\t')
	sample_col = findfirst(headers .== "NOTES") + 1
	samples = headers[sample_col:end]; S = length(samples);

	# The matrix stores genotypes for every sample:
	# 0 = unknown, 1 = homz ref, 2 = hetz, 3 = homz alt
	println("Constructing genotype matrix...")
	genotypes = zeros(Int8, 1_000_000, S); r = 0
	for line in eachline(hetz_snp_vcf)
		c = split(rstrip(line, '\n'), '\t')
		r += 1
		for (s, cell) in enumerate(c[sample_col:end])
			parts = split(cell, ':')
			alt = parse(Int, parts[1]); total = parse(Int, parts[2])
			if total < 30; genotypes[r, s] = 0; continue; end
			b = Binomial(total, 0.5)
			if 0.05 <= cdf(b, alt) <= 0.95
				genotypes[r, s] = 2
			elseif alt / total >= 0.8
				genotypes[r, s] = 3
			elseif alt / total <= 0.2
				genotypes[r, s] = 1
			end
		end
	end
	genotypes = genotypes[1:r, :]

	for i in 1:S
		any_matches = false
		for j in 1:S
			if j == i; continue; end
			matched = 0; total = 0;
			for r in 1:size(genotypes, 1)
				if genotypes[r, i] == 0 || genotypes[r, j] == 0; continue; end
				total += 1
				matched += (genotypes[r, i] == genotypes[r, j])
			end
			if matched / total >= 0.95 && total >= 50
				any_matches = true
				if j > i      # Don't print pairs twice
					println("$(samples[i]) <-> $(samples[j]): $matched / $total matched")
				end
			end
		end

		if any_matches == false
			println("$(samples[i]) does not match with any other samples.")
		end
	end
end

function mutation_rate(vcf_file::IO, coverage_histogram_dir; alt_reads=8)
	vcf = read_vcf(vcf_file); S = length(vcf.sample);
	human_chr = vcat(map(x -> "chr$(x)", 1:22), ["chrX", "chrY"]);
	println("SAMPLE\tPROTEIN ALTERING\tSILENT\tPROTEIN ALTERING (RATE)\tSILENT (RATE)")
	for s in 1:S
		info("Analyzing $(vcf.sample[s])...")
		protein_altering = 0; silent = 0;
		protein_altering_rate = 0; silent_rate = 0;
		histogram = zeros(5000)   # Bins are 0, 1, 2, ..., 4999
		genome_len = 0
		histogram_file = "$(coverage_histogram_dir)/$(vcf.sample[s]).tsv"
		if !isfile(histogram_file)
			warn("No coverage histogram found for $(vcf.sample[s])...")
			continue
		end
		d = readdlm(histogram_file, '\t')
		for chr in human_chr
			genome_len += d[findfirst(x -> x == chr, d[:, 1]), 4]
			for r in 1:size(d, 1)
				if d[r, 1] != chr || d[r, 2] >= 5000; continue; end
				histogram[d[r, 2] + 1] += d[r, 3]
			end
		end

		for r in 1:size(vcf, 1)
			if !vcf.star[r, s]; continue; end
			min_depth = ceil(Int, alt_reads / (vcf.alt[r, s] / vcf.total[r, s]))
			amenable_mbs = (genome_len - sum(histogram[1:min_depth-1])) / 1e6
			info("$(vcf.chromosome[r]):$(vcf.position[r]) - amenable $(amenable_mbs) Mb")
			if is_protein_altering(vcf.effect[r])
				protein_altering += 1
				protein_altering_rate += 1 / amenable_mbs
			else
				silent += 1
				silent_rate += 1 / amenable_mbs
			end
		end

		@printf("%s\t%d\t%d\t%.2f\t%.2f\n", vcf.sample[s], protein_altering, silent, protein_altering_rate, silent_rate)
	end
end

function discard_sketchy_silent(vcf_file::IO)
	line = skip_to_header(vcf_file); println(line)
	notes_col = findfirst(split(line, '\t') .== "NOTES")
	effect_col = findfirst(split(line, '\t') .== "EFFECT")
	for line in eachline(vcf_file)
		cols = split(line, '\t')
		if !is_protein_altering(cols[effect_col])
			if contains(cols[notes_col], "INDEL within"); continue; end
			if contains(cols[notes_col], "Mappability < "); continue; end
		end
		println(line)
	end
end

function sort_samples_by_name(vcf_file::IO)
	line = skip_to_header(vcf_file)
	headers = split(line, '\t')
	notes_col = findfirst(headers .== "NOTES")
	samples = headers[notes_col+1:end]
	order = vcat(1:notes_col, notes_col .+ sortperm(samples))
	println(join(headers[order], '\t'))
	for line in eachline(vcf_file)
		cols = split(line, '\t')
		println(join(cols[order], '\t'))
	end
end

subcommands(somatic, germline, above_background, heterozygous_snps,
	keep_samples, discard_samples,
	discard_shallow, discard_if_frequency_above, discard_blacklisted,
	discard_indels, discard_contingent, discard_sketchy_silent,
	mappability, nearby_indels,
	inside, protein_altering, fractions, statistics,
	predict_effect, annotate, build_annotation_database, error_rate,
	group_by_snp_profile, mutation_rate, sort_by_position, sort_samples_by_name)
