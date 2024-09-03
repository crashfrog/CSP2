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
	nextflow run CSP2.nf --out Test_Output/Contamination_Screen --runmode screen --ref_fasta assets/Screen/Assembly/Lab_Control.fasta --fasta assets/Screen/Assembly/Week_42_Assembly.fasta --reads assets/Screen/Reads --forward _1.fq.gz --reverse _2.fq.gz --readext fq.gz

test: config test_screen ## 10 soil samples - check the relatedness of the cultured isolates
	diff -bur Test_Output/Contamination_Screen/ assets/Screen/Output/Contamination_Screen/