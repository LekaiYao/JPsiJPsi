.PHONY: all DY DPhi Minv SY SPt Tpl clean

DY_LIST := delta_y-0-0.5 delta_y-0.5-1 delta_y-1-1.5 delta_y-1.5-2 delta_y-2-2.5 delta_y-2.5-4
DPhi_LIST := delta_phi-0-0.3927 delta_phi-0.3927-0.7854 delta_phi-0.7854-1.1781 delta_phi-1.1781-1.5708 delta_phi-1.5708-1.9635 delta_phi-1.9635-2.3562 delta_phi-2.3562-2.7489 delta_phi-2.7489-3.1416
Minv_LIST := evt_mass-7.5-17.5 evt_mass-17.5-27.5 evt_mass-27.5-37.5 evt_mass-37.5-47.5 evt_mass-47.5-57.5 evt_mass-57.5-67.5 evt_mass-67.5-107.5
SY_LIST := evt_y-0-0.4 evt_y-0.4-0.8 evt_y-0.8-1.2 evt_y-1.2-1.6 evt_y-1.6-2.0 
SPt_LIST := evt_pt-0-5 evt_pt-5-10 evt_pt-10-15 evt_pt-15-20 evt_pt-20-25 evt_pt-25-30 evt_pt-30-35 evt_pt-35-40 evt_pt-40-80

Model_4D_tot.root: WeightData.root
	@root -l -q Fit_4D_tot.cpp > tmp

$(DY_LIST) $(DPhi_LIST) $(Minv_LIST) $(SY_LIST) $(SPt_LIST): Model_4D_tot.root
	$(eval PARTS := $(subst -, ,$@))
	$(eval VAR := $(word 1,$(PARTS)))
	$(eval MIN := $(word 2,$(PARTS)))
	$(eval MAX := $(word 3,$(PARTS)))
	@root -l -q 'Fit_4D_diff.cpp("$(VAR)", $(MIN), $(MAX))' > tmp
	@root -l -q 'Fit_1D.cpp("$(VAR)", $(MIN), $(MAX))' > tmp
	@echo $(VAR): $(MIN)~$(MAX) >> log.txt
	@root -l -q 'Fit_Check.cpp(true)' | tail -1 | tr -d '\n' >> log.txt
	@root -l -q 'Fit_Check.cpp(false)' | tail -1 | tr -d '\n' >> log.txt
	@root -l -q 'Fit_Check.cpp(false)' | tail -1 >> log.txt

DY: $(DY_LIST)
	@echo >> log.txt

DPhi: $(DPhi_LIST)
	@echo >> log.txt

Minv: $(Minv_LIST)
	@echo >> log.txt

SY: $(SY_LIST)
	@echo >> log.txt

SPt: $(SPt_LIST)
	@echo >> log.txt

Tpl:
	@root -l -q Tpl_Fit.cpp | tail -1

all: DY DPhi Minv SY SPt

clean:
	@rm -f tmp
	@echo > log.txt