#include "/home/luigi/Development/ML4NP/Development/ExportSimData/output_to_csv.C"

void convert(){
	TSystemDirectory *workdir = new TSystemDirectory(".",gSystem->WorkingDirectory());
	TList *RootList = workdir->GetListOfFiles();
	for(auto file : *RootList){
		TString filename(file->GetName());
		if(!filename.BeginsWith("output"))
			continue;
		output_to_csv(filename);
	}

}
