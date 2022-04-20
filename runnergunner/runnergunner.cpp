// SalmonMerge.cpp : Defines the entry point for the application.
//

#include "argparse.h"
#include <filesystem>
#include <vector>
#include <fstream>
#include <memory>
#include <stdlib.h>
#include <set>
#include <string>
#include <string_view>

unsigned int fileSystemMaxFilesOpen = 500;

enum class FileType { Salmon, Tab, Either };

struct DataColumn {
	std::string runname = "none";
	size_t colnum = 0;
};

struct InputFileData {
	std::filesystem::path path = "none";
	FileType filetype = FileType::Salmon;
	std::vector<DataColumn> columns;
	std::shared_ptr<std::ifstream> stream;
	std::shared_ptr<char[]> buffer;
};

enum class RunnerOutput { normal, none, printruns, printgenes };

// C code; should be fast
inline void splitLineOnChar(const std::string& instring, const char& delimiter, std::vector<std::string>& tokens, const size_t& sizehint)
{
	tokens.clear();
	tokens.reserve(sizehint);
	const char* _Beg(&instring[0]), * _End(&instring[instring.size()]);
	for (const char* _Ptr = _Beg; _Ptr < _End; ++_Ptr)
	{
		if (*_Ptr == delimiter)
		{
			tokens.push_back(std::string(_Beg, _Ptr));
			_Beg = 1 + _Ptr;
		}
	}
	tokens.push_back(std::string(_Beg, _End));
}

// C code; should be fast
inline void splitLineOnCharSV(const std::string& instring, const char& delimiter, std::vector<std::string_view>& tokens, const size_t& sizehint)
{
	tokens.reserve(sizehint);
	const char* _Beg(&instring[0]), * _End(&instring[instring.size()]);
	for (const char* _Ptr = _Beg; _Ptr < _End; ++_Ptr)
	{
		if (*_Ptr == delimiter)
		{
			tokens.push_back(std::string(_Beg, _Ptr));
			_Beg = 1 + _Ptr;
		}
	}
	tokens.push_back(std::string(_Beg, _End));
}

inline void splitLineOnTabs(const std::string& line, std::vector<std::string>& tokens, const size_t& sizehint) {
	splitLineOnChar(line, '\t', tokens, sizehint);
}

inline void splitLineOnTabsSVT(const std::string& line, std::vector<std::string_view>& tokens, const size_t& sizehint) {
	splitLineOnCharSV(line, '\t', tokens, sizehint);
}

int _checkSalmonFile(InputFileData & file) {
	std::ifstream newFileStream = std::ifstream(file.path);
	if (newFileStream.good()) {
		// Get first line and check that first line matches expectations
		std::string line;
		std::vector<std::string> linesplit;
		if (std::getline(newFileStream, line)) {
			splitLineOnTabs(line, linesplit, 5);
			if (linesplit.size() != 5) { // wrong number of columns
				std::cerr << "File " << file.path << " should have had 5 columns, but actually had " << linesplit.size() << " and is being omitted\n";
				return 0;
			}
			else if (linesplit.at(3) != "TPM") { // TPM not in correct position
				std::cerr << "Third column of file " << file.path << " should have been TPM, but was actually: " << linesplit.at(3) << ". File is being omitted.\n";
				return 0;
			}
			file.columns.clear();
			file.columns.push_back({ linesplit.at(3) , 3});
		}
		else {
			std::cerr << "Could not get first line from file " << file.path << "\n";
			return 0;
		}
	}
	else {
		std::cerr << "File " << file.path << " failed to open.\n";
		return 0;
	}
	return 1;
}

int _checkTabFile(InputFileData & file) {
	std::ifstream newFileStream = std::ifstream(file.path);
	if (newFileStream.good()) {
		// Get first line and check that first line matches expectations
		std::string line;
		std::vector<std::string> linesplit;
		if (std::getline(newFileStream, line)) {
			splitLineOnTabs(line, linesplit, 10);
			const int numcols = linesplit.size();
			if (numcols < 2) { // wrong number of columns
				std::cerr << "File " << file.path << " should have had at least 2 columns, but actually had " << numcols << " and is being omitted\n";
				return 0;
			}
			else if (linesplit.at(0) != "RNA-see TPM data file") { // TPM not in correct position
				std::cerr << "First cell of file " << file.path << " should have been 'RNA-see TPM data file', but was actually: " << linesplit.at(0) << ". File is being omitted.\n";
				return 0;
			}
			file.columns.clear();
			
			for (size_t i = 1; i < numcols; ++i) file.columns.push_back({ linesplit.at(i) , i });	
		}
		else {
			std::cerr << "Could not get first line from file " << file.path << "\n";
			return 0;
		}
	}
	else {
		std::cerr << "File " << file.path << " failed to open.\n";
		return 0;
	}
	return file.columns.size();
}

int _checkFiles(const std::vector<std::filesystem::path> files, std::vector<InputFileData>& invfiles, const FileType filetype = FileType::Either) {
	int runsum = 0;
	for (auto& file : files) {
		InputFileData filedata;
		filedata.path = file;
		int addedruns = 0;
		if ((filetype == FileType::Salmon || filetype == FileType::Either) && (file.extension() == ".sf")) {
			addedruns = _checkSalmonFile(filedata);
			if (addedruns) {
				filedata.filetype = FileType::Salmon;
				invfiles.push_back(filedata);
			}
			else {
				std::cerr << "Invalid Salmon file: " << file << "\n";
			}
		}
		if ((filetype == FileType::Tab || filetype == FileType::Either) && (file.extension() == ".rnatab")) {
			addedruns = _checkTabFile(filedata);
			if (addedruns) {
				filedata.filetype = FileType::Tab;
				invfiles.push_back(filedata);
			}
			else {
				std::cerr << "Invalid RNA-see tab file: " << file << "\n";
			}
		}
		runsum += addedruns;
	}
	return runsum;
}

void _mergeFilesBatch(std::vector<InputFileData>& batch, const std::string& outFilePath, const FileType filetype,
	RunnerOutput specialmode) {

	// Open input files and identify file types
	for (auto& file : batch) {
		const unsigned long bufSize = 1048576; // 1 mb
		std::shared_ptr<char[]> fileBufferPtr(new char[bufSize]); //<char[]>(bufSize);
		std::string name = file.path.stem().string();
		auto newIfstreamPtr = std::make_shared<std::ifstream>(file.path);
		newIfstreamPtr->rdbuf()->pubsetbuf(fileBufferPtr.get(), bufSize);
		if (newIfstreamPtr->good()) {
			file.stream.swap(newIfstreamPtr);
			file.buffer.swap(fileBufferPtr);
		}
		else {
			std::cerr << "File " << file.path << " failed to open.\n";
			std::cerr << "You may be trying to combine more files than your operating system can simultaneously open.\n";
			exit(1);
		}
	}

	// Open and prep output file
	std::ofstream out;
	const unsigned long bufSize = 10485760; // 10 mb
	auto outBuffer = std::make_unique<char[]>(bufSize);
	if (specialmode != RunnerOutput::none) {
		out = std::ofstream(outFilePath, std::ios::trunc);
		out.rdbuf()->pubsetbuf(outBuffer.get(), bufSize);
		if (!out.good()) {
			throw("Failed to open output file for combined Salmon output");
			exit(1);
		}
	}

	// Prepare to merge the files
	std::string fileLine;
	std::vector<std::string> fileLineSplitVec;

	// Print file header line (if not outputting run/gene list
	if (specialmode == RunnerOutput::normal) {
		out << "RNA-see TPM data file"; // Output header line
	}
	
	// Loop through lines
	bool eof = false;
	bool header = true;

	unsigned long lineNo = 1;
	while (!eof) {
		bool firstfileofline = true;
		std::string genename;
		for (auto& file : batch) {

			// Determine if reached end
			if (!std::getline(*(file.stream.get()), fileLine)) {
				if (!firstfileofline) {
					std::cerr << "File " << file.path << " ended prematurely. Aborting combination operation.\n";
					exit(1);
				}
				eof = true;
				break;
			};

			// Split line
			if (file.filetype == FileType::Salmon) splitLineOnTabs(fileLine, fileLineSplitVec, 5);
			else splitLineOnTabs(fileLine, fileLineSplitVec, file.columns.size());

			// If it's not the header, write/check gene names
			if (!header) {
				if (firstfileofline) { // If it's the first file of the line write gene names in first column
					genename = fileLineSplitVec.at(0);
					if ((specialmode != RunnerOutput::printruns) && (specialmode != RunnerOutput::none)) out << genename;
				}
				else if (genename != fileLineSplitVec.at(0)) { // If it's not the first file, check that gene names match at least
					std::cerr << "Gene name mismatch in file " << file.path << ". Expected gene " << genename << " but read gene " << fileLineSplitVec.at(0) << "\n";
					exit(1);
				}
			}
	
			// copy non-gene data
			
			if (file.filetype == FileType::Salmon) {
				if ((specialmode != RunnerOutput::printgenes) && (specialmode != RunnerOutput::none)) {
					if (header) {
						if (specialmode == RunnerOutput::printruns) {
							out << file.path.stem() << '\n'; // For Salmon files, write file name as run names in header
						}
						else {
							out << '\t' << file.path.stem(); // For Salmon files, write file name as run names in header
						}
					}
					else {
						if (specialmode != RunnerOutput::printruns) {
							out << '\t' << fileLineSplitVec.at(3); // If not header, copy data
						}
					}
				}
			}
			else if (file.filetype == FileType::Tab) { // For tab files, copy run names
				if (!file.columns.size() || fileLineSplitVec.size() < (--file.columns.end())->colnum) {
					std::cerr << "File " << file.path << " ended prematurely. Aborting combination operation.\n";
					exit(1);
				}
				for (auto& col : file.columns) {
					if ((specialmode != RunnerOutput::printgenes) && (specialmode != RunnerOutput::none)) {
						if (specialmode == RunnerOutput::printruns) {
							if (header) {
								out << fileLineSplitVec.at(col.colnum) << '\n';
							}
						}
						else {
							out << '\t' << fileLineSplitVec.at(col.colnum);
						}
					}
				}
			
			}
			firstfileofline = false; // Declare no longer first file of line
		}
		if (!header || (specialmode == RunnerOutput::normal)) {
			if (!eof) out << '\n'; // Terminate line
		}
		header = false; // Declare no longer the header
		if (!(lineNo % 1000)) std::cout << "\rProcessed gene " << lineNo << ".";
		lineNo++;
	}

	// Clean up
	for (auto& file : batch) {
		file.stream->close();
	}
	if (specialmode != RunnerOutput::none) {
		out.close();
	}
}

// Merges the specified .tab or .sf files, assuming that .sf files are named after runs
void _mergeFiles(std::vector<InputFileData>& infiles, const std::string& outfile, bool overwrite, const FileType filetype,
	RunnerOutput specialmode) {
	const size_t numfiles = infiles.size();
	if (numfiles < 1) {
		std::cerr << "Insufficient good files to combine.\n";
		exit(1);
	}

	size_t maxCombine = fileSystemMaxFilesOpen * fileSystemMaxFilesOpen;
	if (numfiles > maxCombine) {
		std::cerr << "Trying to combine too many files (can combine " << maxCombine << " files but tried to combine " << infiles.size() << ").\n";
		exit(1);
	}

	if (filetype == FileType::Either) {
		std::cout << "Processing Salmon and RNA-see tab input files into RNA-see tab output file.\n";
	}
	else if (filetype == FileType::Salmon) {
		std::cout << "Processing Salmon input files into RNA-see tab output file.\n";
	}
	else if (filetype == FileType::Tab) {
		std::cout << "Processing Salmon and RNA-see tab input files into RNA-see tab output file.\n";
	}

	if (std::filesystem::exists(outfile) && !overwrite) {
		std::cerr << "Output file already exists\n";
		exit(1);
	}

	// Check if you have duplicate file names
	std::set<std::filesystem::path> filesAdded;
	for (auto& file : infiles) {
		std::filesystem::path nextPath(file.path);
		if (filesAdded.count(nextPath)) {
			std::cerr << "Trying to merge multiple copies of the same input file";
			exit(1);
		}
		filesAdded.insert(nextPath);
	}

	// divide into merge batches if required
	int batches = ((int)numfiles / fileSystemMaxFilesOpen) + 1;

	if (batches > 1) {
		std::vector<InputFileData> batchtempfiles;
		for (int i = 0; i < batches; ++i) {
			auto start_range = infiles.begin() + i * fileSystemMaxFilesOpen;
			auto end_range = infiles.end() - 1;
			if (i < (batches - 1)) {
				end_range = infiles.begin() + ((i + 1) * fileSystemMaxFilesOpen) - 1;
			}
			std::vector<InputFileData> filebatch;
			while (start_range++ <= end_range) {
				filebatch.push_back(*start_range);
			}
			std::string batchfile = outfile + "_temp_batch" + std::to_string(i);
			batchtempfiles.push_back(InputFileData());
			batchtempfiles.back().path.assign(batchfile);
			batchtempfiles.back().filetype = FileType::Tab;
			_mergeFiles(filebatch, batchfile, overwrite, filetype, specialmode);
		}
		_mergeFilesBatch(batchtempfiles, outfile, filetype, specialmode);
	}
	else {
		_mergeFilesBatch(infiles, outfile, filetype, specialmode);
	}
}

int _removeRuns(std::vector<InputFileData>& files, const std::vector<std::string>& removalvec, bool removedups) {
	std::set<std::string> removals;
	for (auto& run : removalvec) {
		if (removals.count(run) == 0) {
			removals.insert(run);
		}
	}

	int runsum = 0;
	std::vector<InputFileData> newfiles;
	for (auto& file : files) {
		std::vector<DataColumn> newcolumns;
		for (auto& col : file.columns) {

			// Matches name on remove list?
			auto& name = col.runname;
			if (!removals.count(name)) {
				newcolumns.push_back(col); // keep run
				++runsum;
				if (removedups) removals.insert(name); // If removing dups, add to remove list
			}	
		}
		file.columns = newcolumns;
		if (file.columns.size()) newfiles.push_back(file);  // Still has runs, add to new list
	}
	files = newfiles;
	return runsum;
}

// Merges all .tab or .sf files in a directory, assuming that .sf files are named after runs
void mergeFiles(const std::string& outfile, const std::vector<std::filesystem::path>& files, bool overwrite = false, const FileType filetype = FileType::Either,
	RunnerOutput specialmode = RunnerOutput::normal, std::vector<std::string> & removals = std::vector<std::string>(), bool removedups = false) 
{
	std::vector<InputFileData> goodFiles;
	int runsum = _checkFiles(files, goodFiles, filetype);

	if (removedups || removals.size()) {
		std::cout << "Pre-run removal, was going to merge " << runsum << " runs from " << goodFiles.size() << " files, including:\n";
		for (int i = 0; (i < 3) && (i < goodFiles.size()); ++i) {
			std::cout << "\t" << goodFiles.at(i).path << "\n";
		}
		runsum = _removeRuns(goodFiles, removals, removedups);
		std::cout << "Post-run removal, merging " << runsum << " runs from " << goodFiles.size() << " files, including:\n";
		for (int i = 0; (i < 3) && (i < goodFiles.size()); ++i) {
			std::cout << "\t" << goodFiles.at(i).path << "\n";
		}
	}
	else {
		std::cout << "Merging " << runsum << " runs from " << goodFiles.size() << " files, including:\n";
		for (int i = 0; (i < 3) && (i < goodFiles.size()); ++i) {
			std::cout << "\t" << goodFiles.at(i).path << "\n";
		}
	}

	_mergeFiles(goodFiles, outfile, overwrite, FileType::Either, specialmode);
}

// Gathers and merges all .tab or .sf files in a directory, assuming that .sf files are named after runs
void gatherFiles(const std::string& outfile, const std::filesystem::path& dir = "", bool overwrite = false, const FileType filetype = FileType::Either,
	RunnerOutput specialmode = RunnerOutput::normal, std::vector<std::string>& removals = std::vector<std::string>(), bool removedups = false) 
{
	std::cout << "Gathering and checking files from : " << dir << "\n";
	// Loop over directory contents, making a list of good files
	std::vector<std::filesystem::path> checkFiles;
	std::vector<InputFileData> goodFiles;
	for (auto& file : std::filesystem::directory_iterator(dir)) checkFiles.push_back(file.path().string());
	std::cout << "Directory holds " << checkFiles.size() << " files, including:\n";
	for (int i = 0; (i < 3) && (i < checkFiles.size()); ++i) {
		std::cout << "\t" << checkFiles.at(i) << "\n";
	}
	mergeFiles(outfile, checkFiles, overwrite, filetype, specialmode, removals, removedups);
}

// In backend function:
	// Check that all files have same number of rows = genes, and report number at end
	// Check that all output rows have same number of runs, and report number at end

//// Merge together two tab files
//void mergeTabFiles(const std::string& fileA, const std::string& fileB, const std::string& outFile) {
//	std::cout << "Combining RNA-see tab input files:\n\t" << fileA << "\n\t" << fileB << "\n";
//	mergeFiles(outFile, { fileA, fileB }, FileType::Tab);
//}

void print_help(const argparse::ArgumentParser& parser) {
	std::cout << "\nRNA-see runnergunner\n";
	std::cout << "Copyright(c) 2022- Eric Fedosejevs <eric.fedosejevs@gmail.com>\n\n";
	std::cout << "RNA-see runnergunner performs utility operations on quantified RNA-seq data files as output by Salmon, Kallisto etc.\n";
	std::cout << "It is mainly used to gather and merge these runs for analysis by RNA-see.\n\n";
	std::cout << parser;
}

int main(int argc, char* argv[]) {
	bool fakeargs = false;
	if (fakeargs) {
		std::vector<std::string> args = { "-o", "combinedtest.txt", "--remove", "DRR072887", "--duplicates", "--overwrite" }; //,  , "--nooutput" // " - duplicates", 
		std::vector<char*> nargv;
		nargv.push_back(argv[0]);
		for (std::string& s : args) nargv.push_back(&s[0]);
		nargv.push_back(NULL);
		argc = (int)(nargv.size()) - 1;
		argv = nargv.data();
	}

	// print all command line arguments
	std::cout << argv[0] << " ";
	for (auto i = 1; i < argc; ++i) {
		std::cout << argv[i] << " ";
	}
	std::cout << "\n";

	argparse::ArgumentParser program("runnergunner", "0.1", argparse::default_arguments::version);

	program.add_argument("-h", "--help")
		.action([=](const std::string& s) {
		print_help(program);
			})
		.default_value(false)
		.help("shows help message")
		.implicit_value(true)
		.nargs(0);

	program.add_argument("-o", "--output")
		.required()
		.help("specify the output file");

	program.add_argument("-i", "--input")
		.help("specify individual input files, one per -i flag (otherwise gathers all files in directory)");

	program.add_argument("-x", "--remove")
		.help("remove the specified runs from input files");

	program.add_argument("-p", "--duplicates")
		.default_value(false)
		.implicit_value(true)
		.nargs(0)
		.help("remove runs with duplicate names from input files");

	program.add_argument("-r", "--runs")
		.default_value(false)
		.implicit_value(true)
		.nargs(0)
		.help("checks files, but then instead of merging outputs a list of runs from the specified file(s)");

	program.add_argument("-g", "--genes")
		.default_value(false)
		.implicit_value(true)
		.nargs(0)
		.help("checks files, but then instead of merging, outputs a list of genes from the specified file(s)");

	program.add_argument("-n", "--nooutput")
		.default_value(false)
		.implicit_value(true)
		.nargs(0)
		.help("checks file, and goes through a dry merge run without producing any output");

	program.add_argument("-d", "--dir")
		.default_value(std::filesystem::current_path().string())
		.required()
		.nargs(1)
		.help("directory of files being combined");

	program.add_argument("-t", "--type")
		.default_value(std::string("any"))
		.required()
		.nargs(1)
		.help("restrict accepted input file types (salmon (*.sf), rna-see (*.rnatab), any)");

	program.add_argument("-w", "--overwrite")
		.default_value(false)
		.implicit_value(true)
		.nargs(0)
		.help("overwrite existing output file");

	if (argc > 1) {

		try {
			program.parse_args(argc, argv);
		}
		catch (const std::runtime_error& err) {
			std::cerr << err.what() << std::endl;
			std::cerr << program;
			std::exit(1);
		}

		auto dir = program.get<std::string>("--dir");
		if (dir.back() != '/' || dir.back() != '\\') dir.push_back('/'); // add terminating slash to dir
		auto output = program.get<std::string>("--output");

		auto typestr = program.get<std::string>("--type");
		FileType type = FileType::Either;
		if (typestr == "salmon") {
			type = FileType::Salmon;
		}
		else if (typestr == "rna-see") {
			type = FileType::Tab;
		}
		else if (typestr == "default" || typestr == "any") {
			type = FileType::Either;
		}
		else {
			std::cerr << "Invalid input file type specified: " << typestr << "\n.";
			exit(1);
		}

		RunnerOutput specialmode = RunnerOutput::normal;
		bool removedups = false;
		bool overwrite = false;
		std::vector<std::string> removals;
		if (program.is_used("--runs")) {
			specialmode = RunnerOutput::printruns;
			if (program.is_used("--genes")) {
				std::cerr << "Cannot print both run and gene list to same file\n";
				exit(1);
			}
		}
		else if (program.is_used("--genes")) {
			specialmode = RunnerOutput::printgenes;
		}

		if (program.is_used("--nooutput")) {
			specialmode = RunnerOutput::none;
		}
		
		if (program.is_used("--duplicates")) {
			removedups = true;
		}

		if (program.is_used("--overwrite")) {
			overwrite = true;
		}

		if (program.is_used("--remove")) {
			removals = program.get<std::vector<std::string>>("--remove");
			if (!removals.size()) { // if provided removal files
				std::cerr << "Did not supply runs to remove\n";
				exit(1);
			}
		}

		if (program.is_used("--input")) {
			auto inputs = program.get<std::vector<std::string>>("--input");  // {"a.txt", "b.txt", "c.txt"}
			if (inputs.size()) { // if provided input files, do not gather
				std::vector<std::filesystem::path> fullpaths;
				fullpaths.reserve(inputs.size());
				for (std::string& filename : inputs) {
					if (program.is_used("--dir")) { // if provided input files and dir, append file names to dir, else copy
						fullpaths.push_back(std::filesystem::path(dir + filename)); // dir should already have terminating slash
					}
					else {
						fullpaths.push_back(std::filesystem::path(filename)); // dir should already have terminating slash
					}
				}
				mergeFiles(output, fullpaths, overwrite, type, specialmode, removals, removedups);
			}
		}
		else {
			gatherFiles(output, dir, overwrite, type, specialmode, removals, removedups);
		}
		return 0;
	}

	print_help(program);
	return 0;
}
