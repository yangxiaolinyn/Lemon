LIB_DIR = $(PWD)/lib
SRC_DIR = $(PWD)/src
YNOGK_DIR = $(PWD)/ynogk
#Sample_DIR = $(PWD)/src/sampling
#image_DIR = $(PWD)/src/image
main_DIR = $(PWD)
LEMON_DIR = $(PWD)/lemon
#OUT_DIR = $(PWD)/Release
OUT_DIR = Release
MOD_DIR = module

export SRC_DIR LIB_DIR YNOGK_DIR OUT_DIR LEMON_DIR MOD_DIR

FC = mpiifort 
export FC

bin = exe
sample = sampler
im = imager
export bin sample im
default: $(bin)

$(bin): $(SRC_DIR)/*.*90 
	cd $(SRC_DIR) && make $(bin)

#$(sample): $(SRC_DIR)/*.*90 $(YNOGK_DIR)/*.*90 $(Sample_DIR)/*.*90
#	cd $(Sample_DIR) && $(MAKE) -C $(Sample_DIR) $@

#$(im): $(SRC_DIR)/*.*90 $(YNOGK_DIR)/*.*90 $(image_DIR)/*.*90
#	cd $(image_DIR) && $(MAKE) -C $(image_DIR) $@

#.PHONY: clean

clean:
	rm $(SRC_DIR)/$(MOD_DIR)/*.mod
	rm $(SRC_DIR)/$(OUT_DIR)/*.o
	rm $(main_DIR)/exe
 




