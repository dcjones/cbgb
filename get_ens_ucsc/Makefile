
# A simple makefile for quckly gettitg the annotation files I want.

# The following four variables must be defined in the included file:
# 1. name
# 2. ens_genome_name
# 3. ucsc_genom_name
# 4. ens_version

include config.make

name_lc=${shell echo ${name} | tr A-Z a-z}
ens_url=ftp://ftp.ensembl.org/pub/
ens_to_ucsc=~/prj/cbgb/ensemblChromToUCSC/ensemblChromToUCSC


ens_ucsc_gtf=${name}.${ucsc_genome_name}.${ens_version}.gtf
ens_gtf=${name}.${ens_genome_name}.${ens_version}.gtf

${ens_ucsc_gtf} :
	wget ${ens_url}/release-${ens_version}/gtf/${name_lc}/${ens_gtf}.gz -O - \
		| gzip -d \
		| ${ens_to_ucsc} \
		> $@

all : ${ens_ucscs_gtf}


