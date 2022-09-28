# Making all depends on all the .pdf files being generated
all:
	cd ..; ./convert_jupyter_to_python_and_run.sh Tutorial-ETK_thorn-BaikalETK.ipynb && mv Baikal BaikalVacuum BaikalETK/

.PHONY : clean

clean:
	rm -rf Baikal BaikalVacuum
