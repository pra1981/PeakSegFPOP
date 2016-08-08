PeakSegFPOP.out: PeakSegFPOP three.bedGraph
	./PeakSegFPOP three.bedGraph 0.1 | tee PeakSegFPOP.out
PeakSegFPOP: *.cpp *.h compile
	./compile
