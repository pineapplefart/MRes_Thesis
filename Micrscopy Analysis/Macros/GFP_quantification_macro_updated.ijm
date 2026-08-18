// ==== GFP Quantification Macro (.czi, multi-channel) ====
// Expected folder structure:
//   root/1.1/<condition>/*.czi   <- no GFP (used to set threshold)
//   root/2.1/<condition>/*.czi  <- eGFP
//   root/3.1/<condition>/*.czi  <- codon-optimized GFP
// where <condition> = agar, glass, glass_pbs, lc_t0, lc_t3, lc_t5, lc_t7,
//                      lc_t9.5, lc_t24, plug, tapelift
//
// Each .czi file is multi-channel; GFP is channel C2. This macro opens each
// file with Bio-Formats, splits channels, keeps only C2, and runs threshold
// calculation / measurement on that channel alone.
//
// For each condition, this macro:
//   1. Pools GFP-channel pixel intensities across all 1.1 images in that
//      condition folder
//   2. Computes threshold = mean + N*SD (background from the no-GFP strain)
//   3. Applies that threshold to the GFP channel of every image in the
//      matching 2.1 and 3.1 condition folders
//   4. Runs Analyze Particles (thresholded) on each image and saves a
//      per-strain, per-condition summary CSV (Count, Total Area, %Area, Mean)
//
// OUTPUT: one CSV per strain+condition combo, e.g. "2.1_agar_summary.csv",
// saved directly in the root folder you select.
//
// Edit the variables below before running if needed.

thresholdMultiplier = 2;   // threshold = mean + N * SD of the 1.1 (no-GFP) images
refStrain = "1.1";
testStrains = newArray("2.1", "3.1");
gfpChannelPrefix = "C2";   // GFP is channel 2 in these .czi files
fileExtension = ".czi";

rootDir = getDirectory("Choose the root experiment folder (contains 1.1, 2.1, 3.1)");

run("Set Measurements...", "area mean min limit display redirect=None decimal=3");

// ---- Helper: open a .czi file via Bio-Formats and isolate the GFP channel ----
// Leaves only the GFP-channel image open and active; closes all other channels.
function openGFPChannel(filePath, channelPrefix) {
    run("Bio-Formats Importer", "open=[" + filePath + "] color_mode=Default view=Hyperstack stack_order=XYCZT");
    getDimensions(w, h, channels, slices, frames);

    if (channels > 1) {
        run("Split Channels");
        titles = getList("image.titles");
        gfpTitle = "";
        for (i = 0; i < titles.length; i++) {
            if (startsWith(titles[i], channelPrefix + "-")) {
                gfpTitle = titles[i];
            } else {
                selectWindow(titles[i]);
                close();
            }
        }
        if (gfpTitle != "") {
            selectWindow(gfpTitle);
        }
    }
    // if channels == 1, the single opened image is already active -- nothing to do
}

condList = getFileList(rootDir + refStrain + File.separator);

for (c = 0; c < condList.length; c++) {
    condition = condList[c];
    if (!File.isDirectory(rootDir + refStrain + File.separator + condition)) continue;
    condition = replace(condition, "/", "");
    condition = replace(condition, "\\", "");

    refDir = rootDir + refStrain + File.separator + condition + File.separator;
    refFiles = getFileList(refDir);

    sum = 0; sumSq = 0; nPix = 0; imgCount = 0;

    for (f = 0; f < refFiles.length; f++) {
        if (!endsWith(refFiles[f], fileExtension)) continue;
        openGFPChannel(refDir + refFiles[f], gfpChannelPrefix);
        getStatistics(area, mean, min, max, std);
        w = getWidth(); h = getHeight();
        pix = w * h;
        sum += mean * pix;
        sumSq += (std*std + mean*mean) * pix;
        nPix += pix;
        imgCount++;
        close();
    }

    if (imgCount == 0) {
        print("No reference (1.1) images found for condition: " + condition + " -- skipping.");
        continue;
    }

    grandMean = sum / nPix;
    grandVar = (sumSq / nPix) - (grandMean * grandMean);
    grandSD = sqrt(grandVar);
    thresholdVal = grandMean + thresholdMultiplier * grandSD;

    print("Condition [" + condition + "]  n=" + imgCount + " images  mean=" + grandMean + "  SD=" + grandSD + "  threshold=" + thresholdVal);

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
            if (!endsWith(testFiles[f], fileExtension)) continue;
            openGFPChannel(testDir + testFiles[f], gfpChannelPrefix);
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
