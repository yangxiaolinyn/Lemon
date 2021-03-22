LIB_DIR = $(PWD)/lib
SRC_DIR = $(PWD)/src
YNOGK_DIR = $(PWD)/ynogk
LEMON_DIR = $(PWD)/lemon
#OUT_DIR = $(PWD)/Release
OUT_DIR = Release
MOD_DIR = module
#Sample_DIR = $(PWD)/src/sampling
#image_DIR = $(PWD)/src/image
main_DIR = $(PWD)

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
 

#.PHONY: clean

clean:
	rm $(SRC_DIR)/$(MOD_DIR)/*.mod
	rm $(SRC_DIR)/$(OUT_DIR)/*.o
	rm $(main_DIR)/exe 
 




