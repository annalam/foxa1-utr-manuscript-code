
# ADAPTER TRIMMING AND ALIGNMENT
echo *_1.fq.gz | parallel -n1 '[ ! -e ../alignments/${x/_1.fq.gz/.bam} ] &&
mkfifo ${x/.fq.gz/.reads.pipe} && 
mkfifo ${x/_1.fq.gz/_2.reads.pipe} && 
(cutadapt -f fastq -m 20 -a AGATCGGAAGAGC -A AGATCGGAAGAGC $x ${x/_1.fq/_2.fq} -o >(fasta trim by quality - 30 | fasta mask by quality - 20 | fasta to raw - > ${x/.fq.gz/.reads.pipe}) -p >(fasta trim by quality - 30 | fasta mask by quality - 20 | fasta to raw - > ${x/_1.fq.gz/_2.reads.pipe}) &) && 
bowtie2 -r -p8 -X 1000 --score-min L,0,-0.6 -x ~/tools/bowtie2-indexes/homo_sapiens/hg38 -1 ${x/.fq.gz/.reads.pipe} -2 ${x/_1.fq.gz/_2.reads.pipe} | samblaster | samtools view -u - | samtools sort -@8 -o ../alignments/${x/_1.fq.gz/.bam} -; 
rm ${x/.fq.gz/.reads.pipe} ${x/_1.fq.gz/_2.reads.pipe}'
echo *.bam | parallel -n8 '[ ! -e ${x}.bai ] && samtools index $x'







# ESTIMATE PANEL SIZE BY CALCULATING COVERAGE HISTOGRAM FOR ENTIRE GENOME
echo *-UMI.bam | parallel -n8 'samtools view -b -q10 -F 0x400 $x | bedtools genomecov -ibam stdin -g ~/homo_sapiens/hg38.chrom.sizes > ../coverage_histograms/${x/.bam/.tsv}'

tsv_files = readdir();
histogram = zeros(2000, length(tsv_files));
total = zeros(1, length(tsv_files));
human_chr = vcat(map(x -> "chr$(x)", 1:22), ["chrX", "chrY"]);
for (s, tsv_file) in enumerate(tsv_files)
	d = readtsv(tsv_file)
	for chr in human_chr
		total[s] += d[findfirst(x -> x == chr, d[:, 1]), 4]
		for r in find((d[:, 1] .== chr) & (d[:, 2] .< 2000))
			histogram[d[r, 2] + 1, s] += d[r, 3]
		end
	end
end

threshold = 100;
for s in 1:length(tsv_files)
	captured = round(Int, total[s] - sum(histogram[1:threshold, s]))
	println("$(tsv_files[s]): $(captured)")
end
median(map(s -> total[s] - sum(histogram[1:threshold, s]),
	1:length(tsv_files)))
total[1] - sum(median(histogram, 2)[1:threshold])









# IDENTIFY SOMATIC MUTATIONS AND GERMLINE VARIANTS IN TARGETED PANEL
echo X `seq 22` Y | parallel -n5 'mutato call2 --region=chr${x} --alt-reads=8 --alt-frac=0.005 --min-mapq=0 ~/homo_sapiens/hg38.fa *.bam > ../mutations/chr${x}.vcf'
cat chr1.vcf <(cat chr{2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y}.vcf | grep -v 'CHROM\|##') > variants.vcf
variant nearby indels variants.vcf | variant somatic --alt-reads=10 --alt-frac=0.01 --test-ref-ratio=3 --test-bg-ratio=25 --ref-reads=20 --min-sidedness=15 --min-mapq=10 - ../tumor_normal_pairs.txt | variant predict effect - | variant mappability - ~/homo_sapiens/mappability_170bp_hg38.jld 90 | variant discard sketchy silent - | variant annotate - ~/homo_sapiens/cosmic_77_hg38.jld | variant annotate - ~/homo_sapiens/exac_0.3_hg38.jld | variant discard if frequency above - ~/homo_sapiens/exac_0.3_hg38.jld 0.005 > somatic.vcf






# SUPPLEMENTARY TABLE 1: LIQUID BIOPSIES
using Distributions
d = readtsv("table_s4_mutations.tsv");
headers = d[1, :][:]; d = d[2:end, :];
samples = unique(d[:, 1]);
chromosome = d[:, findone(headers .== "Chrom")];
total_reads = d[:, findone(headers .== "Total reads")];
maf = [isempty(x) ? NaN : parse(Float64, x[1:end-1]) / 100
	for x in d[:, findone(headers .== "Allele-%")]];
effect = d[:, findone(headers .== "Effect")];

# Estimate ctDNA fractions using code borrowed from "clonality" tool
cancer_frac = zeros(length(samples));
for s in 1:length(samples)
	for r in find(d[:, 1] .== samples[s])
		if isnan(maf[r]); continue; end

		ploidy = 2
		if endswith(chromosome[r], 'X') || endswith(chromosome[r], 'Y')
			ploidy = 1
		end

		# The observed mutant allele fraction is drawn from a binomial
		# distribution Bin(n, p) where p = true mutant allele fraction.
		# Here we see what the true MAF could be, assuming that the
		# observed MAF is a 95% quantile anomaly.
		true_maf = NaN;
		for p in 0:0.005:1
			dist = Binomial(total_reads[r], p)
			if cdf(dist, round(Int, maf[r] * total_reads[r])) < 0.95
				true_maf = p; break;
			end
		end

		if ploidy == 1
			this_cancer_frac = true_maf
		elseif ploidy == 2
			# Since we cannot detect copy number changes when cancer-%
			# is low, we conservatively always assume LOH for autosomes.
			this_cancer_frac = 2 / (1 / true_maf + 1)
		end

		cancer_frac[s] = max(cancer_frac[s], this_cancer_frac)
	end
end

# Output as file "ctdna_fractions.tsv"
out = open("ctdna_fractions.tsv", "w");
for s in 1:length(samples)
	@printf(out, "%s\t%.3f\n", samples[s], cancer_frac[s])
end
close(out);

# Count somatic mutations of different types
total_somatic = zeros(length(samples));
total_somatic_cds = zeros(length(samples));
total_somatic_5utr = zeros(length(samples));
total_somatic_3utr = zeros(length(samples));
for s in 1:length(samples)
	for r in find(d[:, 1] .== samples[s])
		total_somatic[s] += 1
		if is_coding_region(effect[r])
			total_somatic_cds[s] += 1
		elseif rx(effect[r], "5'-UTR")
			total_somatic_5utr[s] += 1
		elseif rx(effect[r], "3'-UTR")
			total_somatic_3utr[s] += 1
		end
	end
end

# Perform some sanity checks on the table
for s in 1:length(samples)
	if cancer_frac[s] == 0 && total_somatic[s] > 0
		info("$(samples[s]) has rearrangements but no mutations.")
	end
end

# How many total in coding regions, UTRs, and intronic/intergenic?
sum(total_somatic_cds)
sum(total_somatic_5utr) + sum(total_somatic_3utr)
2000 - sum(total_somatic_cds) - sum(total_somatic_5utr) - sum(total_somatic_3utr)








# RESULTS: CALCULATE SOMATIC MUTATION RATES

# Figure 1B barplot
with_ctdna = 439;
total_cds_len = sum([r.len for r in regions if r.class == "CDS"]) / 1e6;
total_utr5_len = sum(covered_len[region_type .== "5'-UTR"]) / 1e6;
total_utr3_len = sum(covered_len[region_type .== "3'-UTR"]) / 1e6;

cds_mutations / total_cds_len / with_ctdna
utr5_mutations / total_utr5_len / with_ctdna
utr3_mutations / total_utr3_len / with_ctdna
(utr5_mutations + utr3_mutations) / (total_utr5_len + total_utr3_len) / with_ctdna

using Distributions
1 - cdf(Binomial(cds_mutations + utr5_mutations + utr3_mutations, total_cds_len / (total_cds_len + total_utr5_len + total_utr3_len)), cds_mutations)

# What fraction substitutions, what fraction indels?
count(m -> m.region == "CDS" && m.indel, mutations) / count(m -> m.region == "CDS", mutations)
count(m -> endswith(m.region, "UTR") && m.indel, mutations) / count(m -> endswith(m.region, "UTR"), mutations)

using Statistics
valid = [m.region == "CDS" || endswith(m.region, "UTR") for m in mutations];
fisher_exact([m.region == "CDS" for m in mutations[valid]],
	[m.indel for m in mutations[valid]])









# FIGURE 1C AND 1D: BARPLOTS OF MUTATION COUNT
cd ~/datasets/foxa1_utr
using Plot2;
with_ctdna = 439;

# Figure 1C
d = readtsv("table_s4_mutations.tsv"); d = d[2:end, :];
cds_truncating = [fill("", 0) for g in 1:G];
cds_missense = [fill("", 0) for g in 1:G];
utr5_mutations = [fill("", 0) for g in 1:G];
utr3_mutations = [fill("", 0) for g in 1:G];
for r in 1:size(d, 1)
	sample = d[r, 1]
	g = in_region(d[r, :], "5'-UTR")
	if g > 0; push!(utr5_mutations[g], sample); end

	g = in_region(d[r, :], "3'-UTR")
	if g > 0; push!(utr3_mutations[g], sample); end

	g = findone(genes .== d[r, 6])
	if g > 0
		if is_truncating(d[r, 7])
			push!(cds_truncating[g], sample)
		elseif is_missense(d[r, 7])
			push!(cds_missense[g], sample)
		end
	end
end

order = sortperm(length.(cds_truncating) + length.(cds_missense) + length.(utr5_mutations) + length.(utr3_mutations), rev=true)[1:15];
stacked_bar_plot(hcat(length.(cds_truncating), length.(cds_missense), length.(utr5_mutations), length.(utr3_mutations))[order, :])

for g in order
	@printf("%s\t%.1f%%\n", genes[g], length(unique(vcat(cds_truncating[g], cds_missense[g], utr5_mutations[g], utr3_mutations[g]))) / with_ctdna * 100)
end

# Figure 1D
with_ctdna = 439;
genes = unique([r.gene for r in regions]);
utr3_len = [sum(r.covered_len for r in regions if r.gene == g && r.class == "3'-UTR") for g in genes];
utr3_subs = [count(m -> m.gene == g && m.in_utr3 && !m.indel, mutations)
	for g in genes];
utr3_indels = [count(m -> m.gene == g && m.in_utr3 && m.indel, mutations)
	for g in genes];

utr3_subs_rate = utr3_subs ./ (utr3_len ./ 1e6) ./ with_ctdna;
utr3_indels_rate = utr3_indels ./ (utr3_len ./ 1e6) ./ with_ctdna;

order = sortperm(utr3_subs_rate + utr3_indels_rate, rev=true)[1:15];
using Plot2; stacked_bar_plot(hcat(utr3_subs_rate, utr3_indels_rate)[order, :])

for g in order
	@printf("%s\t%.1f%%\n", genes[g], length(unique(vcat(utr3_subs[g], utr3_indels[g]))) / with_ctdna * 100)
end






# FIGURE 1E: MUTATION LOLLIPOP PLOT
cd ~/datasets/foxa1_utr
# Mutations from our own cohort
foxa1_muts = [m for m in mutations if m.gene == "FOXA1" && m.primary];
position = [round(Int, m.position + m.mut_len / 2) for m in foxa1_muts];
ref_allele = [m.ref_allele for m in foxa1_muts];
alt_allele = [m.alt_allele for m in foxa1_muts];

# Mutations from cBioPortal
d = readtsv("cbioportal_mutations.tsv"); d = d[2:end, :];
d = d[d[:, 2] .== "chr14", :];
position = Float64.(d[:, 3]);
ref_allele = d[:, 4]; alt_allele = d[:, 5];
for k in 1:length(ref_allele)
	if ref_allele[k] == "-"
		ref_allele[k] = "N"; alt_allele[k] = "N$(alt_allele[k])";
	elseif alt_allele[k] == "-"
		ref_allele[k] = "N$(ref_allele[k])"; alt_allele[k] = "N";
	end
end

# Calculate the dot positions in the scatter plot
using Plot2
mut_mrna_pos = zeros(Int, 0); mut_color = Vector{RGB}(0);
for r in 1:length(position)
	# Calculate mutation's position inside the mRNA transcript
	mrna_pos = 0
	if 37_594_901 <= position[r] <= 37_595_120         # First exon (220 bp)
		mrna_pos = 37_595_120 - position[r] + 1
	elseif 37_589_552 <= position[r] <= 37_592_711     # Second exon (3160 bp)
		mrna_pos = 37_592_711 - position[r] + 220 + 1
	else
		continue
	end

	push!(mut_mrna_pos, mrna_pos)
	if is_indel(ref_allele[r], alt_allele[r])
		if mod(length(alt_allele[r]) - length(ref_allele[r]), 3) == 0
			if length(alt_allele[r]) > length(ref_allele[r])
				push!(mut_color, RGB(255, 128, 128))
			else
				push!(mut_color, RGB(128, 128, 255))
			end
		else
			if length(alt_allele[r]) > length(ref_allele[r])
				push!(mut_color, RGB(255, 0, 0))
			else
				push!(mut_color, RGB(0, 0, 255))
			end
		end
	else
		push!(mut_color, RGB(0, 255, 0))
	end
end

order = sortperm(mut_mrna_pos);
mut_mrna_pos = mut_mrna_pos[order];
mut_color = mut_color[order];

# Stagger the mutations on the y-axis to make them not overlap
mrna_len = 3160 + 220;
y_level = fill(NaN, length(mut_mrna_pos));
y_space = falses(50, mrna_len);
for k in 1:length(mut_mrna_pos)
	level = findfirst(y_space[:, mut_mrna_pos[k]] .== false)
	window = max(mut_mrna_pos[k] - 20, 1):min(mut_mrna_pos[k] + 20, mrna_len)
	y_space[level, window] = true
	y_level[k] = level
end

figure("~/plot.pdf", size=(10, 3)) do
	scatter_plot(mut_mrna_pos, y_level, color=mut_color, marker="o", size=11)
	xlim(0, 3380); ylim(0, 40);
end







# FIGURE 2B: BASE COMPOSITION AROUND FOXA1 3'-UTR INDELS
using BioSequences, Variant
genome_path = expanduser("~/homo_sapiens/hg38.fa");
genome = FASTA.Reader(open(genome_path), index="$(genome_path).fai");
chr14 = sequence(genome["chr14"]);
foxa1_utr3_start = 37589551; foxa1_utr3_end = 37591364;
d = readtsv("table_s4_mutations.tsv"); d = d[2:end, :];
junctions = fill("", 0);
for r in 1:size(d, 1)
	if d[r, 6] != "FOXA1"; continue; end
	if !is_indel(d[r, 4], d[r, 5]); continue; end
	if any(k -> d[k, 6] == "FOXA1" && d[k, 3] == d[r, 3], 1:r-1); continue; end
	pos = Int(d[r, 3])
	if foxa1_utr3_start <= pos <= foxa1_utr3_end
		junction = "$(chr14[pos-19:pos])|$(chr14[pos+1:pos+20])"
		push!(junctions, junction)
		#println("$(junction)\t$(d[r, 4])\t$(d[r, 5])")
	end
end

gc_frac = fill(NaN, 20);
for k in 1:20
	gc_frac[k] = count(j -> j[21 - k] in "GC", junctions) + count(j -> j[21 + k] in "GC", junctions)
end
gc_frac /= 2 * length(junctions)
using Plot2; line_plot(1:length(gc_frac), gc_frac)

# Generate list of random 10 bp sequences from within the 3'-UTR
junctions = fill("", 0);
for k in 1:1000
	pos = rand(foxa1_utr3_start+50:foxa1_utr3_end-50)
	push!(junctions, "$(chr14[pos-19:pos])|$(chr14[pos+1:pos+20])")
end








# FIGURE 2C: SEQUENCE MOTIFS AROUND FOXA1 3'-UTR INDELS
using BioSequences, Variant, Plot2
genome_path = expanduser("~/homo_sapiens/hg38.fa");
genome = FASTA.Reader(open(genome_path), index="$(genome_path).fai");
chr14 = ReferenceSequence(sequence(genome["chr14"]));
foxa1_utr3_start = 37589551; foxa1_utr3_end = 37591364;

K = 4;
kmers = [DNAKmer{K}(UInt64(k)) for k in 0:(4^K-1)];

# Regions near somatic indels
kmer_count_near_indels = zeros(length(kmers)); J = 0;
d = readtsv("table_s4_mutations.tsv"); d = d[2:end, :];
for r in 1:size(d, 1)
	if d[r, 6] != "FOXA1"; continue; end
	if !is_indel(d[r, 4], d[r, 5]); continue; end
	if any(k -> d[k, 6] == "FOXA1" && d[k, 3] == d[r, 3], 1:r-1); continue; end
	pos = Int(d[r, 3])
	if !(foxa1_utr3_start <= pos <= foxa1_utr3_end); continue; end
	junction = chr14[pos-19:pos+20]
	for (_, kmer) in each(DNAKmer{K}, junction)
		kmer_count_near_indels[convert(UInt64, canonical(kmer)) + 1] += 1
		J += 1
	end
end
kmer_count_near_indels /= J;

# All regions inside FOXA1 3'-UTR
kmer_count_everywhere = zeros(length(kmers)); J = 0;
for pos in foxa1_utr3_start:(foxa1_utr3_end - 39)
	junction = chr14[pos:pos+39]
	for (_, kmer) in each(DNAKmer{K}, junction)
		kmer_count_everywhere[convert(UInt64, canonical(kmer)) + 1] += 1
		J += 1
	end
end
kmer_count_everywhere /= J;

figure("~/plot.pdf", size=(3.5,3)) do
	valid = (kmer_count_everywhere .> 0) | (kmer_count_near_indels .> 0)
	scatter_plot(kmer_count_everywhere[valid], kmer_count_near_indels[valid], marker="o", point_size=20)
	if K == 3; scale = 0.12; xlim(0, scale); ylim(0, scale); end
	if K == 4; scale = 0.035; xlim(0, scale); ylim(0, scale); end
end


order = sortperm(kmer_count_near_indels, rev=true);
hcat(kmers, kmer_count_near_indels, kmer_count_everywhere)[order[1:15], :]








# ARE FOXA1 MUTATIONS SUBCLONAL OR CLONAL?
d = readtsv("ctdna_fractions.tsv");
ctdna_frac = Dict{String, Float64}(d[r, 1] => d[r, 2] for r in 1:size(d, 1));
# Run first lines of code from Figure 1C code here...

foxa1_cds_mut = zeros(0); foxa1_utr3_mut = zeros(0);
other_cds_mut = zeros(0); other_utr3_mut = zeros(0);
for m in [m for m in mutations if m.chromosome != "chrX"]
	ctdna = ctdna_frac[m.sample]
	caf = m.alt_frac / ctdna
	if m.in_cds
		push!(m.gene == "FOXA1" ? foxa1_cds_mut : other_cds_mut, caf)
	end
	if m.in_utr3
		push!(m.gene == "FOXA1" ? foxa1_utr3_mut : other_utr3_mut, caf)
	end
end

using Plot2, Statistics;
figure("~/plot.pdf", size=(3,3)) do
	beeswarm_plot(foxa1_cds_mut, other_cds_mut, foxa1_utr3_mut, other_utr3_mut, point_size=8)
end
ranksum_test(foxa1_cds_mut, other_cds_mut)
ranksum_test(foxa1_utr3_mut, other_utr3_mut)

ranksum_test(vcat(foxa1_cds_mut, other_cds_mut), vcat(foxa1_utr3_mut, other_utr3_mut))








# SUPPLEMENTARY FIGURE: CO-OCCURRENCE OF FOXA1 AND OTHER RECURRENT MUTATIONS
using Variant
d = readtsv("table_s4_mutations.tsv"); d = d[2:end, :];
samples = unique(d[:, 1]); S = length(samples);

foxa1_cds_mutation = falses(S);
foxa1_utr3_mutation = falses(S);
tp53_mutation = falses(S);
ctnnb1_mutation = falses(S);
apc_mutation = falses(S);
pik3ca_mutation = falses(S);
pten_mutation = falses(S);
akt1_mutation = falses(S);
spop_mutation = falses(S);
ar_mutation = falses(S);

for r in 1:size(d, 1)
	gene = d[r, 6]; effect = d[r, 7];
	s = findone(samples .== d[r, 1])
	if gene == "FOXA1" && is_coding_region(effect)
		foxa1_cds_mutation[s] = true
	elseif gene == "FOXA1" && contains(effect, "3'-UTR")
		foxa1_utr3_mutation[s] = true
	elseif gene == "TP53" && is_protein_altering(effect)
		tp53_mutation[s] = true
	elseif gene == "CTNNB1" && is_protein_altering(effect)
		ctnnb1_mutation[s] = true
	elseif gene == "APC" && is_truncating(effect)
		apc_mutation[s] = true
	elseif gene == "PIK3CA" && is_protein_altering(effect)
		pik3ca_mutation[s] = true
	elseif gene == "PTEN" && is_protein_altering(effect)
		pten_mutation[s] = true
	elseif gene == "AKT1" && is_protein_altering(effect)
		akt1_mutation[s] = true
	elseif gene == "SPOP" && is_protein_altering(effect)
		spop_mutation[s] = true
	elseif gene == "AR" && is_protein_altering(effect) && d[r, 3] >= 67711620
		ar_mutation[s] = true
	end
end

using Plot2
order = hierarchical_order(foxa1_utr3_mutation, foxa1_cds_mutation, tp53_mutation, ctnnb1_mutation, apc_mutation, pik3ca_mutation, pten_mutation, akt1_mutation, spop_mutation, ar_mutation);
matrix = MatrixLayer(10, S, glyph=GLYPH_TILE);
shade(x) = x ? RGB(0, 0, 0) : RGB(230, 230, 230);
matrix.color[1, :] = map(shade, foxa1_utr3_mutation);
matrix.color[2, :] = map(shade, foxa1_cds_mutation);
matrix.color[3, :] = map(shade, tp53_mutation);
matrix.color[4, :] = map(shade, ctnnb1_mutation);
matrix.color[5, :] = map(shade, apc_mutation);
matrix.color[6, :] = map(shade, pik3ca_mutation);
matrix.color[7, :] = map(shade, pten_mutation);
matrix.color[8, :] = map(shade, akt1_mutation);
matrix.color[9, :] = map(shade, spop_mutation);
matrix.color[10, :] = map(shade, ar_mutation);
matrix_plot([matrix[:, order]], path="~/plot.svg", cell_width=26)

using Statistics
@printf("%.3f\n", fisher_exact(foxa1_utr3_mutation, foxa1_cds_mutation));
@printf("%.3f\n", fisher_exact(foxa1_utr3_mutation, tp53_mutation));
@printf("%.3f\n", fisher_exact(foxa1_utr3_mutation, ctnnb1_mutation));
@printf("%.3f\n", fisher_exact(foxa1_utr3_mutation, apc_mutation));
@printf("%.3f\n", fisher_exact(foxa1_utr3_mutation, pik3ca_mutation));
@printf("%.3f\n", fisher_exact(foxa1_utr3_mutation, pten_mutation));
@printf("%.3f\n", fisher_exact(foxa1_utr3_mutation, akt1_mutation));
@printf("%.3f\n", fisher_exact(foxa1_utr3_mutation, spop_mutation));
@printf("%.3f\n", fisher_exact(foxa1_utr3_mutation, ar_mutation));

ref = foxa1_utr3_mutation | foxa1_cds_mutation;
@printf("%.3f\n", fisher_exact(ref, tp53_mutation));
@printf("%.3f\n", fisher_exact(ref, ctnnb1_mutation));
@printf("%.3f\n", fisher_exact(ref, apc_mutation));
@printf("%.3f\n", fisher_exact(ref, pik3ca_mutation));
@printf("%.3f\n", fisher_exact(ref, pten_mutation));
@printf("%.3f\n", fisher_exact(ref, akt1_mutation));
@printf("%.3f\n", fisher_exact(ref, spop_mutation));
@printf("%.3f\n", fisher_exact(ref, ar_mutation));

sum(foxa1_utr3_mutation | foxa1_cds_mutation)
sum(foxa1_cds_mutation)




