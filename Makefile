BuildDir := build

all: $(BuildDir) | copl_build

copl_build: $(BuildDir)
	@make --no-print-directory -C src
	@make --no-print-directory -C rpt
	@echo run build/copl

clean:
	@rm -rf $(BuildDir)
	@echo clean up build
	
$(BuildDir):
	@mkdir $@
