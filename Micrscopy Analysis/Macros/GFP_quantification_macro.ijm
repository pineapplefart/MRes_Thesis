// ==== GFP Quantification Macro ====
// Expected folder structure:
//   root/1.1/<condition>/*.tif   <- no GFP (used to set threshold)
//   root/2.1/<condition>/*.tif   <- eGFP
//   root/3.1/<condition>/*.tif   <- codon-optimized GFP
// where <condition> = agar, glass, glass_pbs, lc_t0, lc_t3, lc_t5, lc_t7,
//                      lc_t9.5, lc_t24, plug, tapelift
//
// For each condition, this macro:
//   1. Pools pixel intensities across all 1.1 images in that condition folder
//   2. Computes threshold = mean + N*SD (background from the no-GFP strain)
//   3. Applies that threshold to every image in the matching 2.1 and 3.1
//      condition folders
//   4. Runs Analyze Particles (thresholded) on each image and saves a
//      per-strain, per-condition summary CSV (Count, Total Area, %Area, Mean)
//
// OUTPUT: one CSV per strain+condition combo, e.g. "2.1_agar_summary.csv",
// saved directly in the root folder you select.
//
// Edit the two variables below before running if needed.

thresholdMultiplier = 3; // threshold = mean + N * SD of the 1.1 (no-GFP) images
refStrain = "1.1";
testStrains = newArray("2.1", "3.1");

rootDir = getDirectory("Choose the root experiment folder (contains 1.1, 2.1, 3.1)");

run("Set Measurements...", "area mean min limit display redirect=None decimal=3");

condList = getFileList(rootDir + refStrain + File.separator);

for (c = 0; c < condList.length; c++) {
    condition = condList[c];
    if (!File.isDirectory(rootDir + refStrain + File.separator + condition)) continue;
    condition = replace(condition, "/", "");
    condition = replace(condition, "\\", "");

    refDir = rootDir + refStrain + File.separator + condition + File.separator;
    refFiles = getFileList(refDir);

    sum = 0; sumSq = 0; nPix = 0; nImages = 0;

    for (f = 0; f < refFiles.length; f++) {
        if (!(endsWith(refFiles[f], ".tif") || endsWith(refFiles[f], ".tiff"))) continue;
        open(refDir + refFiles[f]);
        getStatistics(area, mean, min, max, std);
        w = getWidth(); h = getHeight();
        pix = w * h;
        sum += mean * pix;
        sumSq += (std*std + mean*mean) * pix;
        nPix += pix;
        nImages++;
        close();
    }

    if (nImages == 0) {
        print("No reference (1.1) images found for condition: " + condition + " -- skipping.");
        continue;
    }

    grandMean = sum / nPix;
    grandVar = (sumSq / nPix) - (grandMean * grandMean);
    grandSD = sqrt(grandVar);
    thresholdVal = grandMean + thresholdMultiplier * grandSD;

    print("Condition [" + condition + "]  n=" + nImages + " images  mean=" + grandMean + "  SD=" + grandSD + "  threshold=" + thresholdVal);

    for (s = 0; s < testStrains.length; s++) {
        strain = testStrains[s];
        testDir = rootDir + strain + File.separator + condition + File.separator;
        if (!File.exists(testDir)) {
            print("  Missing matching folder for strain " + strain + ": " + testDir);
            continue;
        }
        testFiles = getFileList(testDir);

        run("Clear Results");
        if (isOpen("Summary")) {
            selectWindow("Summary");
            run("Close");
        }

        nTestImages = 0;
        for (f = 0; f < testFiles.length; f++) {
            if (!(endsWith(testFiles[f], ".tif") || endsWith(testFiles[f], ".tiff"))) continue;
            open(testDir + testFiles[f]);
            setThreshold(thresholdVal, 1e9);
            run("Analyze Particles...", "size=0-Infinity display clear summarize");
            close();
            nTestImages++;
        }

        if (nTestImages == 0) {
            print("  No images found in: " + testDir);
            continue;
        }

        if (isOpen("Summary")) {
            selectWindow("Summary");
            saveAs("Results", rootDir + strain + "_" + condition + "_summary.csv");
            run("Close");
        }
    }
}

print("Done. Per-strain, per-condition summary CSVs saved in: " + rootDir);
