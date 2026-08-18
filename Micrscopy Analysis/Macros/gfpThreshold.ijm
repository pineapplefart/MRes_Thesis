setThreshold(2547, 65535);
setOption("BlackBackground", true);
run("Convert to Mask");
run("Analyze Particles...", "size=0-Infinity display summarize");