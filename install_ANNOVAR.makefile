
# ~~~~~ ANNOVAR ~~~~~ # 
# RUN cd /opt/bin && wget http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.revision150617.tar.gz && tar xvfz annovar.revision150617.tar.gz && rm -f annovar.revision150617.tar.gz
# RUN mkdir -p "${ANNOVAR_DB_DIR}"
# RUN "${ANNOVAR_DIR}/annotate_variation.pl" -downdb -buildver hg19 -webfrom annovar clinvar_20150629 "${ANNOVAR_DB_DIR}"
# RUN "${ANNOVAR_DIR}/annotate_variation.pl" -downdb -buildver hg19 -webfrom annovar cytoBand "${ANNOVAR_DB_DIR}"
ANNOVAR_DIR:=./bin/annovar
ANNOVAR_DB_DIR:=./bin/annovar/db/hg19

annovar.revision150617.tar.gz:
	wget http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.revision150617.tar.gz

$(ANNOVAR_DB_DIR)/hg19_ALL.sites.2015_08.txt:
	"${ANNOVAR_DIR}/annotate_variation.pl" -downdb -buildver hg19 -webfrom annovar 1000g2015aug "${ANNOVAR_DB_DIR}"

1000g2015aug: $(ANNOVAR_DB_DIR)/hg19_ALL.sites.2015_08.txt

$(ANNOVAR_DB_DIR)/hg19_refGene.txt:
	"${ANNOVAR_DIR}/annotate_variation.pl" -downdb -buildver hg19 -webfrom annovar refGene "${ANNOVAR_DB_DIR}" 

refGene: $(ANNOVAR_DB_DIR)/hg19_refGene.txt

$(ANNOVAR_DB_DIR)/hg19_cosmic70.txt:
	"${ANNOVAR_DIR}/annotate_variation.pl" -downdb -buildver hg19 -webfrom annovar cosmic70 "${ANNOVAR_DB_DIR}"

cosmic70: $(ANNOVAR_DB_DIR)/hg19_cosmic70.txt

$(ANNOVAR_DB_DIR)/hg19_clinvar_20170905.txt:
	"${ANNOVAR_DIR}/annotate_variation.pl" -downdb -buildver hg19 -webfrom annovar clinvar_20170905 "${ANNOVAR_DB_DIR}"

clinvar_20170905: $(ANNOVAR_DB_DIR)/hg19_clinvar_20170905.txt

$(ANNOVAR_DB_DIR)/hg19_intervar_20170202.txt:
	"${ANNOVAR_DIR}/annotate_variation.pl" -downdb -buildver hg19 -webfrom annovar intervar_20170202 "${ANNOVAR_DB_DIR}"

intervar_20170202: $(ANNOVAR_DB_DIR)/hg19_intervar_20170202.txt


$(ANNOVAR_DB_DIR)/hg19_dbnsfp33a.txt:
	"${ANNOVAR_DIR}/annotate_variation.pl" -downdb -buildver hg19 -webfrom annovar dbnsfp33a "${ANNOVAR_DB_DIR}"
	
dbnsfp33a: $(ANNOVAR_DB_DIR)/hg19_dbnsfp33a.txt

$(ANNOVAR_DB_DIR)/hg19_esp6500siv2_all.txt:
	"${ANNOVAR_DIR}/annotate_variation.pl" -downdb -buildver hg19 -webfrom annovar esp6500siv2_all "${ANNOVAR_DB_DIR}"

esp6500siv2_all: $(ANNOVAR_DB_DIR)/hg19_esp6500siv2_all.txt

$(ANNOVAR_DB_DIR)/hg19_kaviar_20150923.txt:
	"${ANNOVAR_DIR}/annotate_variation.pl" -downdb -buildver hg19 -webfrom annovar kaviar_20150923 "${ANNOVAR_DB_DIR}"

kaviar_20150923: $(ANNOVAR_DB_DIR)/hg19_kaviar_20150923.txt


$(ANNOVAR_DB_DIR)/hg19_gnomad_exome.txt:
	"${ANNOVAR_DIR}/annotate_variation.pl" -downdb -buildver hg19 -webfrom annovar gnomad_exome "${ANNOVAR_DB_DIR}"

gnomad_exome: $(ANNOVAR_DB_DIR)/hg19_gnomad_exome.txt

$(ANNOVAR_DB_DIR)/hg19_gnomad_genome.txt: 
	"${ANNOVAR_DIR}/annotate_variation.pl" -downdb -buildver hg19 -webfrom annovar gnomad_genome "${ANNOVAR_DB_DIR}"

gnomad_genome: $(ANNOVAR_DB_DIR)/hg19_gnomad_genome.txt

$(ANNOVAR_DB_DIR)/hg19_avsnp150.txt:
	"${ANNOVAR_DIR}/annotate_variation.pl" -downdb -buildver hg19 -webfrom annovar avsnp150 "${ANNOVAR_DB_DIR}"

avsnp150: $(ANNOVAR_DB_DIR)/hg19_avsnp150.txt

$(ANNOVAR_DB_DIR)/hg19_cadd13gt10.txt:
	"${ANNOVAR_DIR}/annotate_variation.pl" -downdb -buildver hg19 -webfrom annovar cadd13gt10 "${ANNOVAR_DB_DIR}"

cadd13gt10: $(ANNOVAR_DB_DIR)/hg19_cadd13gt10.txt

$(ANNOVAR_DB_DIR)/hg19_fathmm.txt:
	"${ANNOVAR_DIR}/annotate_variation.pl" -downdb -buildver hg19 -webfrom annovar fathmm "${ANNOVAR_DB_DIR}"

fathmm: $(ANNOVAR_DB_DIR)/hg19_fathmm.txt

$(ANNOVAR_DB_DIR)/hg19_eigen.txt.gz: 
	"${ANNOVAR_DIR}/annotate_variation.pl" -downdb -buildver hg19 -webfrom annovar eigen "${ANNOVAR_DB_DIR}"

eigen: $(ANNOVAR_DB_DIR)/hg19_eigen.txt

annovar_db: 1000g2015aug refGene clinvar_20170905 intervar_20170202 dbnsfp33a esp6500siv2_all kaviar_20150923 gnomad_exome gnomad_genome avsnp150 fathmm eigen

