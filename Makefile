four.out: PeakSegFPOP four.bedGraph
	./PeakSegFPOP four.bedGraph 0.1 | tee four.out
three.out: PeakSegFPOP three.bedGraph
	./PeakSegFPOP three.bedGraph 0.1 | tee three.out
PeakSegFPOP: *.cpp *.h compile
	./compile
