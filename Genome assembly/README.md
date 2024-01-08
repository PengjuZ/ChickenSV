Genome assembly:

The genome size of Wenchang chicken was estimated using Jellyfish and GenomeScope based on k-mers obtained from Illumina paired-end reads. The genome assembly pipeline involved various steps including Hi-C based binning, hybrid assembly, ONT assembly, scaffolding, gap filling, and polishing. Phased contig assemblies were generated using HiFi and Hi-C reads, and hybrid assemblies were assembled using PacBio HiFi and ONT ultra-long reads. The resulting contigs were anchored to the chromosome level using Salsa2 and YaHS, and manual corrections were made using Juicebox assembly tools. The scaffolds were further scaffolded and gap-closed using RagTag and manual methods, and then iteratively polished using HiFi and Illumina paired-end reads.

Genome quality assessment:

The assembly quality was evaluated using the vgp-assembly pipeline and the Asset evaluation tool. The telomeric-identifier tool was used to identify telomere regions. Genome completeness was assessed using the BUSCO program, and base accuracy was measured using Merqury.
