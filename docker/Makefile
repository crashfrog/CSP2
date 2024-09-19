.PHONY: help

.ONESHELL:

help: ## Show this help
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}'

config:
	@ cat << EOF
	profiles {
		standard {
			process.executor = 'local'
			params.cores = `nproc --all`
		}
	}
	EOF > ~/.nextflow/config


ifeq (screen, $(firstword $(MAKECMDGOALS)))
	runargs := $(wordlist 2, $(words $(MAKECMDGOALS)), $(MAKECMDGOALS))
	$(eval $(runargs):;@true)
endif

ifeq (snp, $(firstword $(MAKECMDGOALS)))
	runargs := $(wordlist 2, $(words $(MAKECMDGOALS)), $(MAKECMDGOALS))
	$(eval $(runargs):;@true)
endif

screen: config ## determine whether query isolates are close to a reference

snp: config ## generate pairwise distances for a set of query isolates

snpdiffs: config

test_screen:
	nextflow run CSP2.nf -profile standard --runmode screen --fasta assets/Screen/Assembly/Week_42_Assembly.fasta --reads assets/Screen/Reads/ --ref_fasta assets/Screen/Assembly/Lab_Control.fasta --out ./CSP2_Test_Screen --readext fq.gz --forward _1.fq.gz --reverse _2.fq.gz

test_snp:
	nextflow run CSP2.nf -profile standard --runmode snp --fasta assets/SNP/ --n_ref 3 --out ./CSP2_Test_SNP --max_missing 50

test: config test_screen test_snp ## 10 soil samples - check the relatedness of the cultured isolates
	ls -lah assets/Screen/Output/Contamination_Screen/
	diff -bur ./CSP2_Test_SNP/snpdiffs assets/SNP/Output/Soil_Analysis/snpdiffs