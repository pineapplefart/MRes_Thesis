// ============================================================================
//  Re-analysis from saved ROIs  --  GFP quantification project
// ============================================================================
//  Re-measures images using the ROI sets you already drew, so NO clicking is
//  needed. Runs in batch mode (windows hidden) and is therefore fast.
//
//  WHAT IT DOES
//    Pass 1 -- for each condition, opens the 1.1 (no-GFP) images, restricts to
//              their saved ROIs, and pools the GFP-channel pixel statistics to
//              get  threshold = mean + N x SD   (N is set below; you asked for 2).
//    Pass 2 -- opens the 2.1 / 3.1 images of that condition, loads their saved
//              ROIs, and measures each cell against that condition's threshold.
//
//  Pooling is done properly across images: sums and sums-of-squares are
//  accumulated over every ROI pixel in the condition, so a big image does not
//  count the same as a small one. (ImageJ reports the sample SD, n-1; at these
//  pixel counts the difference from the population SD is negligible.)
//
//  INPUTS
//    - the image root folder (with subfolders on, to walk strain/condition)
//    - the folder holding the previously saved <name>_ROIs.zip files
//    - an output folder for the new CSVs
//
//  OUTPUTS  (suffixed with the multiplier, so 3-SD results are never overwritten)
//    reanalysis_cell_measurements_<N>SD.csv -- one row per cell
//    reanalysis_thresholds_<N>SD.csv        -- the threshold used per condition
//
//  Images with no saved ROI set are skipped and listed in the Log.
// ============================================================================

// ---------------- Parameters (all editable in the startup dialog) -----------
var fileExtension        = ".czi";

// Only accept filenames of the form 202607XX_<strain>_tX_-XX.czi (X = digits).
// The strain group is 2-4 numbers separated by dashes OR underscores, so
// 3-1-1, 1_1_2 and 1_1 all match.
// Full-string match, so anything with an extra token (e.g. _1500ms_) is excluded.
// The final number may be written with or without a leading dash (_-01 or _01),
// and \p{Pd} accepts any dash character, so en dashes pasted from elsewhere
// still match. The timepoint may be a whole number (t3), a decimal (t9.5), or a
// decimal split across tokens (t9_5). To allow months other than July, change
// 202607 to 2026\d{2}.
var filenameRegex        = "202607\\d{2}_\\d+(?:[\\p{Pd}_]\\d+){1,3}_t\\d+(?:[._]\\d+)?_\\p{Pd}?\\d+\\.czi";
var enforceFilenameRegex = true;

var conditionFilter      = "lc_";      // only conditions starting with this
var useConditionFilter   = true;
var excludeConditions    = "lc_t24";   // comma-separated condition folders to drop
var excludeTimepoints    = "24";       // comma-separated t-values to drop by filename token

var thresholdMultiplier  = 2;          // threshold = mean + N x SD   (was 3)

var backgroundStrain     = "1.1";      // no-GFP reference strain
var targetStrains        = "2.1,3.1";  // strains to measure

var gfpChannelPrefix     = "C2";
var brightfieldPrefix    = "C1";
var measureInPixels      = true;
var listFilesAndStop     = false;      // dry run: list what passes the filters, then quit

var importerCommand      = "Bio-Formats Windowless Importer";

// ---------------- Startup dialog --------------------------------------------
Dialog.create("Re-analyse from saved ROIs");
Dialog.addMessage("Reuses your existing ROI sets. No clicking required.");
Dialog.addNumber("Threshold multiplier (mean + N x SD):", thresholdMultiplier);
Dialog.addString("File extension:", fileExtension);
Dialog.addCheckbox("Require filename pattern", enforceFilenameRegex);
Dialog.addString("   pattern (regex):", filenameRegex, 44);
Dialog.addCheckbox("Only conditions starting with:", useConditionFilter);
Dialog.addString("   condition prefix:", conditionFilter);
Dialog.addString("Exclude conditions:", excludeConditions);
Dialog.addString("Exclude timepoints (t values):", excludeTimepoints);
Dialog.addString("Background strain:", backgroundStrain);
Dialog.addString("Target strains:", targetStrains);
Dialog.addCheckbox("Measure in pixels (ignore micron calibration)", measureInPixels);
Dialog.addCheckbox("Dry run: list matched/excluded files, then stop", listFilesAndStop);
Dialog.show();

thresholdMultiplier  = Dialog.getNumber();
fileExtension        = Dialog.getString();
enforceFilenameRegex = Dialog.getCheckbox();
filenameRegex        = Dialog.getString();
useConditionFilter   = Dialog.getCheckbox();
conditionFilter      = Dialog.getString();
excludeConditions    = Dialog.getString();
excludeTimepoints    = Dialog.getString();
backgroundStrain     = Dialog.getString();
targetStrains        = Dialog.getString();
measureInPixels      = Dialog.getCheckbox();
listFilesAndStop     = Dialog.getCheckbox();

// NOTE: these must be declared with "var". In the ImageJ macro language a plain
// top-level assignment is NOT visible inside functions, and roiPathFor() needs
// roiDir. Without var the interpreter treats it as an undefined numeric value.
var imageRoot = "";
var roiDir    = "";
var outputDir = "";

imageRoot = getDirectory("Choose the IMAGE root folder (contains 1.1 / 2.1 / 3.1)");
roiDir    = getDirectory("Choose the folder containing the saved _ROIs.zip files");
outputDir = getDirectory("Choose an output folder for the new CSVs");

nTag       = "" + thresholdMultiplier;
csvPath    = outputDir + "reanalysis_cell_measurements_" + nTag + "SD.csv";
threshPath = outputDir + "reanalysis_thresholds_" + nTag + "SD.csv";

// ---------------- Globals ----------------------------------------------------
var collected     = newArray(0);
var excludedLog   = newArray(0);
var condNames     = newArray(0);   // condition folder names, e.g. lc_t3
var condN         = newArray(0);   // pixel count pooled from 1.1
var condSum       = newArray(0);   // sum of pixel values
var condSumSq     = newArray(0);   // sum of squared pixel values
var condThreshold = newArray(0);
var condImages    = newArray(0);   // how many 1.1 images contributed

var openedC1 = "";
var openedC2 = "";

// ---------------- Helpers ----------------------------------------------------
function stripExtension(name) {
    if (endsWith(toLowerCase(name), toLowerCase(fileExtension)))
        return substring(name, 0, lengthOf(name) - lengthOf(fileExtension));
    return name;
}

function inCsvList(value, csvList) {
    parts = split(csvList, ",");
    for (i = 0; i < parts.length; i++) {
        if (toLowerCase(String.trim(parts[i])) == toLowerCase(String.trim(value))) return true;
    }
    return false;
}

// Pull the t-value out of a filename. Handles "_t3_", "_t9.5_" and "_t9_5_"
// (a decimal split across two tokens). For the split form the following token is
// only joined on if it is all digits AND is not the last token, so the trailing
// image number in a name like "..._t3_01.czi" is never mistaken for a decimal.
function timepointOf(name) {
    parts = split(stripExtension(name), "_");
    for (i = 0; i < parts.length; i++) {
        p = parts[i];
        if (startsWith(toLowerCase(p), "t") && lengthOf(p) > 1) {
            rest = substring(p, 1);
            if (matches(rest, "[0-9.]+")) {
                if (i + 1 < parts.length - 1 && matches(parts[i + 1], "\\d+"))
                    return rest + "." + parts[i + 1];
                return rest;
            }
        }
    }
    return "";
}

function noteExcluded(path, reason) {
    // Include the detected strain/condition so a wrong root folder is obvious:
    // if "strain" is not 1.1/2.1/3.1, the root folder is at the wrong level.
    ctx = File.getName(File.getParent(File.getParent(path))) + "/"
        + File.getName(File.getParent(path)) + "/";
    excludedLog = Array.concat(excludedLog, ctx + File.getName(path) + "   [" + reason + "]");
}

function fileAllowed(path) {
    name = File.getName(path);
    cond = File.getName(File.getParent(path));

    if (useConditionFilter && conditionFilter != "") {
        if (!startsWith(toLowerCase(cond), toLowerCase(conditionFilter))) {
            noteExcluded(path, "condition " + cond + " not " + conditionFilter + "*");
            return false;
        }
    }
    if (excludeConditions != "" && inCsvList(cond, excludeConditions)) {
        noteExcluded(path, "condition " + cond + " excluded");
        return false;
    }
    tp = timepointOf(name);
    if (excludeTimepoints != "" && tp != "" && inCsvList(tp, excludeTimepoints)) {
        noteExcluded(path, "timepoint t" + tp + " excluded");
        return false;
    }
    if (enforceFilenameRegex && !matches(name, filenameRegex)) {
        noteExcluded(path, "filename does not match the pattern");
        return false;
    }
    return true;
}

function collectFiles(dir) {
    list = getFileList(dir);
    for (i = 0; i < list.length; i++) {
        if (endsWith(list[i], "/")) {
            collectFiles(dir + list[i]);
        } else if (endsWith(toLowerCase(list[i]), toLowerCase(fileExtension))) {
            if (fileAllowed(dir + list[i])) collected = Array.concat(collected, dir + list[i]);
        }
    }
}

function strainOf(path)    { return File.getName(File.getParent(File.getParent(path))); }
function conditionOf(path) { return File.getName(File.getParent(path)); }

function roiPathFor(path) {
    return roiDir + stripExtension(File.getName(path)) + "_ROIs.zip";
}

function condIndex(name) {
    for (i = 0; i < condNames.length; i++) { if (condNames[i] == name) return i; }
    return -1;
}

function ensureCond(name) {
    idx = condIndex(name);
    if (idx >= 0) return idx;
    condNames     = Array.concat(condNames, name);
    condN         = Array.concat(condN, 0);
    condSum       = Array.concat(condSum, 0);
    condSumSq     = Array.concat(condSumSq, 0);
    condThreshold = Array.concat(condThreshold, -1);
    condImages    = Array.concat(condImages, 0);
    return condNames.length - 1;
}

function reduceToMiddleSlice(windowTitle) {
    selectWindow(windowTitle);
    total = nSlices();
    if (total > 1) {
        mid = floor((total + 1) / 2);
        run("Duplicate...", "title=" + windowTitle + "_mid duplicate range=" + mid + "-" + mid);
        close(windowTitle);
        rename(windowTitle);
    }
}

function openAndSplitChannels(path) {
    openedC1 = "";
    openedC2 = "";
    before = getList("image.titles");
    run(importerCommand, "open=[" + path + "] autoscale color_mode=Default view=Hyperstack stack_order=XYCZT");
    after = getList("image.titles");

    fresh = newArray(0);
    for (i = 0; i < after.length; i++) {
        wasOpen = false;
        for (j = 0; j < before.length; j++) { if (after[i] == before[j]) wasOpen = true; }
        if (!wasOpen) fresh = Array.concat(fresh, after[i]);
    }
    if (fresh.length == 0) return;

    selectWindow(fresh[0]);
    getDimensions(w0, h0, ch, sl, fr);
    if (fresh.length == 1 && ch > 1) {
        selectWindow(fresh[0]);
        run("Split Channels");
        after = getList("image.titles");
        fresh = newArray(0);
        for (i = 0; i < after.length; i++) {
            wasOpen = false;
            for (j = 0; j < before.length; j++) { if (after[i] == before[j]) wasOpen = true; }
            if (!wasOpen) fresh = Array.concat(fresh, after[i]);
        }
    } else if (fresh.length == 1 && ch == 1) {
        selectWindow(fresh[0]);
        close();
        return;
    }

    leftovers = newArray(0);
    for (i = 0; i < fresh.length; i++) {
        if (startsWith(fresh[i], brightfieldPrefix + "-") || startsWith(fresh[i], brightfieldPrefix + "_")) openedC1 = fresh[i];
        else if (startsWith(fresh[i], gfpChannelPrefix + "-") || startsWith(fresh[i], gfpChannelPrefix + "_")) openedC2 = fresh[i];
        else leftovers = Array.concat(leftovers, fresh[i]);
    }
    if (openedC1 == "" && openedC2 == "" && leftovers.length == 2) {
        openedC1 = leftovers[0];
        openedC2 = leftovers[1];
        leftovers = newArray(0);
        print("   note: no C1-/C2- prefixes; assuming window order (1st = C1, 2nd = C2).");
    }
    for (i = 0; i < leftovers.length; i++) {
        if (isOpen(leftovers[i])) { selectWindow(leftovers[i]); close(); }
    }
    if (openedC1 != "") reduceToMiddleSlice(openedC1);
    if (openedC2 != "") reduceToMiddleSlice(openedC2);
    if (measureInPixels) {
        if (openedC1 != "") { selectWindow(openedC1); run("Set Scale...", "distance=0 known=0 unit=pixel"); }
        if (openedC2 != "") { selectWindow(openedC2); run("Set Scale...", "distance=0 known=0 unit=pixel"); }
    }
}

function upperLimitForBitDepth() {
    bd = bitDepth();
    if (bd == 8)  return 255;
    if (bd == 16) return 65535;
    return 1e30;
}

function closeIfOpen(t) { if (t != "" && isOpen(t)) { selectWindow(t); close(); } }

// ---------------- Setup -------------------------------------------------------
run("Set Measurements...", "area mean standard min median integrated area_fraction display redirect=None decimal=4");
roiManager("reset");
run("Clear Results");
print("\\Clear");
if (measureInPixels) units = "pixels"; else units = "calibrated";

collectFiles(imageRoot);
print("Re-analysis with threshold = mean + " + thresholdMultiplier + " x SD");
print("Matched " + collected.length + " file(s); excluded " + excludedLog.length + ".");

if (listFilesAndStop) {
    print("");
    print("--- MATCHED (" + collected.length + ") ---");
    for (i = 0; i < collected.length; i++) {
        marker = "";
        if (!File.exists(roiPathFor(collected[i]))) marker = "   [NO SAVED ROIs]";
        print("  " + strainOf(collected[i]) + "/" + conditionOf(collected[i])
              + "/" + File.getName(collected[i]) + marker);
    }
    print("");
    print("--- EXCLUDED (" + excludedLog.length + ") ---");
    print("  pattern in use: " + filenameRegex);
    print("  (shown as detected_strain/detected_condition/filename [reason])");
    print("  If the strain column is not " + backgroundStrain + " / " + targetStrains
          + ", you picked the wrong image root folder.");
    for (i = 0; i < excludedLog.length; i++) print("  " + excludedLog[i]);
    exit("Dry run complete -- see the Log. Nothing was measured.");
}

if (collected.length == 0) exit("No files matched the filters. Try the dry run to see why.");

// Write headers only if the files are new, so a re-run after an error does not
// insert a second header row. Note that re-running APPENDS: delete or rename the
// old CSVs first if you want a clean file.
if (!File.exists(threshPath))
    File.append("Condition,NumBackgroundImages,PixelsPooled,BackgroundMean,BackgroundSD,Multiplier,Threshold", threshPath);
else
    print("NOTE: " + threshPath + " already exists -- new rows will be appended.");
if (!File.exists(csvPath))
    File.append("File,Strain,Condition,Timepoint,CellID,SelectionMethod,Area,GFP_Mean,GFP_SD,GFP_Median,GFP_Max,GFP_RawIntDen,GFP_PctAreaAboveThreshold,GFP_Threshold,Multiplier,Units", csvPath);
else
    print("NOTE: " + csvPath + " already exists -- new rows will be appended.");

setBatchMode(true);
missingRois = 0;

// ---------------- PASS 1: pool background stats from the no-GFP strain --------
print("");
print("=== Pass 1: pooling " + backgroundStrain + " background statistics ===");
for (f = 0; f < collected.length; f++) {
    filePath = collected[f];
    // Assign to a variable before comparing. The macro interpreter assumes a
    // user function returns a number when the call appears directly inside a
    // comparison, which raises "Numeric return value expected".
    thisStrain = strainOf(filePath);
    if (thisStrain != backgroundStrain) continue;

    roiPath = roiPathFor(filePath);
    if (!File.exists(roiPath)) {
        print("   no ROIs for " + File.getName(filePath) + " -- not contributing to the threshold");
        missingRois++;
        continue;
    }

    openAndSplitChannels(filePath);
    if (openedC2 == "") { print("   could not open GFP channel: " + File.getName(filePath)); continue; }

    roiManager("reset");
    roiManager("Open", roiPath);
    ci = ensureCond(conditionOf(filePath));

    selectWindow(openedC2);
    for (r = 0; r < roiManager("count"); r++) {
        selectWindow(openedC2);
        roiManager("Select", r);
        nPix = getValue("Area");          // pixels, calibration removed
        m    = getValue("Mean");
        sd   = getValue("StdDev");
        condN[ci]     = condN[ci] + nPix;
        condSum[ci]   = condSum[ci] + m * nPix;
        condSumSq[ci] = condSumSq[ci] + nPix * (sd * sd + m * m);
    }
    condImages[ci] = condImages[ci] + 1;

    roiManager("reset");
    closeIfOpen(openedC1);
    closeIfOpen(openedC2);
}

// ---------------- Compute one threshold per condition -------------------------
print("");
print("=== Thresholds (mean + " + thresholdMultiplier + " x SD) ===");
for (i = 0; i < condNames.length; i++) {
    if (condN[i] <= 1) {
        condThreshold[i] = -1;
        print("   " + condNames[i] + ": no usable background pixels -- images in this condition will be SKIPPED");
        continue;
    }
    bgMean = condSum[i] / condN[i];
    bgVar  = condSumSq[i] / condN[i] - bgMean * bgMean;
    if (bgVar < 0) bgVar = 0;
    bgSD   = sqrt(bgVar);
    condThreshold[i] = bgMean + thresholdMultiplier * bgSD;
    print("   " + condNames[i] + ": mean " + d2s(bgMean, 2) + ", SD " + d2s(bgSD, 2)
          + "  ->  threshold " + d2s(condThreshold[i], 2)
          + "   (" + condImages[i] + " image(s), " + d2s(condN[i], 0) + " px)");
    File.append(condNames[i] + "," + condImages[i] + "," + d2s(condN[i], 0) + ","
        + d2s(bgMean, 4) + "," + d2s(bgSD, 4) + "," + thresholdMultiplier + ","
        + d2s(condThreshold[i], 4), threshPath);
}

// ---------------- PASS 2: measure the target strains --------------------------
print("");
print("=== Pass 2: measuring " + targetStrains + " ===");
cellsWritten = 0;
for (f = 0; f < collected.length; f++) {
    filePath = collected[f];
    strain   = strainOf(filePath);
    if (!inCsvList(strain, targetStrains)) continue;

    cond = conditionOf(filePath);
    ci   = condIndex(cond);
    if (ci < 0 || condThreshold[ci] < 0) {
        print("   skipping " + File.getName(filePath) + " -- no threshold for " + cond);
        continue;
    }
    roiPath = roiPathFor(filePath);
    if (!File.exists(roiPath)) {
        print("   no ROIs for " + File.getName(filePath) + " -- skipped");
        missingRois++;
        continue;
    }

    openAndSplitChannels(filePath);
    if (openedC2 == "") { print("   could not open GFP channel: " + File.getName(filePath)); continue; }

    roiManager("reset");
    roiManager("Open", roiPath);
    thr = condThreshold[ci];
    tp  = timepointOf(File.getName(filePath));

    selectWindow(openedC2);
    upper = upperLimitForBitDepth();
    for (r = 0; r < roiManager("count"); r++) {
        selectWindow(openedC2);
        roiManager("Select", r);
        roiName = Roi.getName;
        if (startsWith(roiName, "manual")) method = "drawn"; else method = "clicked";

        resetThreshold();
        a      = getValue("Area");
        m      = getValue("Mean");
        sd     = getValue("StdDev");
        median = getValue("Median");
        mx     = getValue("Max");
        intden = getValue("RawIntDen");

        setThreshold(thr, upper);
        pctAbove = getValue("%Area");
        resetThreshold();

        File.append(File.getName(filePath) + "," + strain + "," + cond + "," + tp + ","
            + roiName + "," + method + "," + d2s(a, 4) + "," + d2s(m, 4) + ","
            + d2s(sd, 4) + "," + d2s(median, 4) + "," + d2s(mx, 4) + ","
            + d2s(intden, 4) + "," + d2s(pctAbove, 4) + "," + d2s(thr, 4) + ","
            + thresholdMultiplier + "," + units, csvPath);
        cellsWritten++;
    }
    print("   " + strain + "/" + cond + "/" + File.getName(filePath)
          + ": " + roiManager("count") + " cell(s)");

    roiManager("reset");
    closeIfOpen(openedC1);
    closeIfOpen(openedC2);
}

setBatchMode(false);
roiManager("reset");
print("");
print("----------------------------------------------------------");
print("Done. " + cellsWritten + " cell row(s) written.");
if (missingRois > 0) print("WARNING: " + missingRois + " file(s) had no saved ROI set and were skipped.");
print("Per-cell measurements: " + csvPath);
print("Thresholds used:       " + threshPath);
