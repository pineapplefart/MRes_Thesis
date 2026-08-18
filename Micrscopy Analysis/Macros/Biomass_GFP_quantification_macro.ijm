// ==== GFP Quantification Macro (.czi, multi-channel, biomass-normalized) ====
// Expected folder structure:
//   root/1.1/<condition>/*.czi   <- no GFP (used to set threshold)
//   root/2.1/<condition>/*.czi  <- eGFP
//   root/3.1/<condition>/*.czi  <- codon-optimized GFP
// where <condition> = agar, glass, glass_pbs, lc_t0, lc_t3, lc_t5, lc_t7,
//                      lc_t9.5, lc_t24, plug, tapelift
//
// Each .czi file has (at least) 2 channels: C1 = brightfield, C2 = GFP.
//
// For every image (regardless of strain/condition), this macro:
//   1. Segments C1 (brightfield) into a "biomass" mask -- any fungal material
//      present, whether discrete cells/spores or continuous mycelium. No
//      per-cell segmentation is attempted, so this works the same way across
//      conditions where cells are discrete (glass, glass_pbs) and conditions
//      where they form mycelium (agar, plug, tapelift, later timepoints).
//   2. Restricts ALL subsequent GFP (C2) measurements to that biomass region,
//      so images with more or less fungal material in frame are normalized
//      against each other rather than compared as raw whole-image values.
//
// For each condition, the workflow is:
//   a. Pool GFP-channel pixel intensities from WITHIN THE BIOMASS REGION of
//      all 1.1 (no-GFP) images in that condition folder.
//   b. Compute threshold = mean + N*SD of that pooled biomass-restricted
//      background.
//   c. Apply that threshold to the biomass-restricted GFP channel of every
//      2.1 and 3.1 image in the matching condition folder.
//   d. Save a per-strain, per-condition summary CSV with, per image:
//        Count            - number of distinct GFP+ regions within biomass
//        Total Area       - GFP+ pixel area within biomass
//        %Area            - GFP+ area AS A PERCENT OF BIOMASS AREA (not
//                            whole image), since Analyze Particles computes
//                            %Area relative to the active selection
//        Mean             - mean intensity of just the GFP+ pixels
//        BiomassArea      - total biomass pixel area (from C1 mask) for that
//                            image, i.e. the normalization denominator
//        BiomassMeanIntensity - mean GFP intensity across ALL biomass pixels
//                            (not just thresholded ones) -- a continuous,
//                            threshold-independent readout of GFP per unit
//                            biomass
//
// IMPORTANT -- VERIFY BEFORE BATCH RUNNING:
// Brightfield segmentation direction depends on your imaging setup. Test on
// one image first (e.g. run just the makeBiomassSelection portion manually,
// or watch the first file processed) to confirm the mask actually selects
// fungal material and not empty background. If the mask looks inverted,
// flip biomassDarkBackground below.

thresholdMultiplier = 3;
refStrain = "1.1";
testStrains = newArray("2.1", "3.1");
fileExtension = ".czi";

// Brightfield (C1) segmentation settings for biomass masking.
biomassThresholdMethod = "Default";
biomassDarkBackground = false; // set true if fungal material appears BRIGHTER than background in your brightfield images; leave false if it appears DARKER (typical for transmitted light)

rootDir = getDirectory("Choose the root experiment folder (contains 1.1, 2.1, 3.1)");

run("Set Measurements...", "area mean min limit display redirect=None decimal=3");

var c1Title, c2Title;
var biomassAreaResult;

// ---- Opens a .czi file via Bio-Formats, splits channels, keeps only C1 and C2 ----
function openBothChannels(filePath) {
    run("Bio-Formats Importer", "open=[" + filePath + "] color_mode=Default view=Hyperstack stack_order=XYCZT");
    getDimensions(w, h, channels, slices, frames);
    c1Title = "";
    c2Title = "";
    if (channels > 1) {
        run("Split Channels");
        titles = getList("image.titles");
        for (i = 0; i < titles.length; i++) {
            if (startsWith(titles[i], "C1-")) {
                c1Title = titles[i];
            } else if (startsWith(titles[i], "C2-")) {
                c2Title = titles[i];
            } else {
                selectWindow(titles[i]);
                close();
            }
        }
    } else {
        c2Title = getTitle();
    }
}

// ---- Segments C1 into a biomass mask, creates a selection, closes C1 ----
// Leaves the selection loaded (retrievable via "Restore Selection") and sets
// the global biomassAreaResult to its pixel area (0 if nothing was detected).
function makeBiomassSelection() {
    selectWindow(c1Title);
    if (biomassDarkBackground) {
        setAutoThreshold(biomassThresholdMethod + " dark");
    } else {
        setAutoThreshold(biomassThresholdMethod);
    }
    run("Create Selection");
    if (selectionType() == -1) {
        biomassAreaResult = 0;
    } else {
        biomassAreaResult = getValue("Area");
    }
    resetThreshold();
    close();
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

        openBothChannels(refDir + refFiles[f]);
        if (c1Title == "" || c2Title == "") {
            print("  Warning: could not find both channels in " + refFiles[f] + " -- skipping.");
            if (isOpen(c1Title)) { selectWindow(c1Title); close(); }
            if (isOpen(c2Title)) { selectWindow(c2Title); close(); }
            continue;
        }

        makeBiomassSelection(); // segments + closes C1, leaves a selection ready to restore

        if (biomassAreaResult <= 0) {
            print("  No biomass detected in " + refFiles[f] + " -- skipping.");
            selectWindow(c2Title);
            close();
            continue;
        }

        selectWindow(c2Title);
        run("Restore Selection");
        getStatistics(area, mean, min, max, std);
        sum += mean * area;
        sumSq += (std*std + mean*mean) * area;
        nPix += area;
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

    print("Condition [" + condition + "]  n=" + imgCount + " images  biomass-mean=" + grandMean + "  SD=" + grandSD + "  threshold=" + thresholdVal);

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

            openBothChannels(testDir + testFiles[f]);
            if (c1Title == "" || c2Title == "") {
                print("  Warning: could not find both channels in " + testFiles[f] + " -- skipping.");
                if (isOpen(c1Title)) { selectWindow(c1Title); close(); }
                if (isOpen(c2Title)) { selectWindow(c2Title); close(); }
                continue;
            }

            makeBiomassSelection();

            if (biomassAreaResult <= 0) {
                print("  No biomass detected in " + testFiles[f] + " -- skipping.");
                selectWindow(c2Title);
                close();
                continue;
            }

            selectWindow(c2Title);
            run("Restore Selection");
            getStatistics(biomassArea, biomassMeanGFP, biomassMin, biomassMax, biomassStd);

            setThreshold(thresholdVal, 1e9);
            run("Analyze Particles...", "size=0-Infinity display clear summarize");

            row = Table.size("Summary") - 1;
            Table.set("BiomassArea", row, biomassArea, "Summary");
            Table.set("BiomassMeanIntensity", row, biomassMeanGFP, "Summary");
            Table.update("Summary");

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
