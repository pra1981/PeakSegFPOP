PeakSegFPOP.out: PeakSegFPOP three.bedGraph
	./PeakSegFPOP three.bedGraph 0.1 | tee PeakSegFPOP.out
PeakSegFPOP: *.cpp *.h
	g++ -std=c++11 funPieceListLog.cpp PeakSegFPOPLog.cpp -o PeakSegFPOP
