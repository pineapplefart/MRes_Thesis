// ============================================================================
//  Interactive Cell Segmentation (click-to-select)  --  GFP quantification project
// ============================================================================
//  Replaces the fully automatic brightfield segmentation for datasets where the
//  cell types are too different for one auto-threshold to handle.
//
//  WORKFLOW PER IMAGE:
//    1. The .czi is opened, split into C1 (brightfield) and C2 (GFP), and
//       reduced to its middle Z-slice if it is a stack.
//    2. Two windows appear: a slightly blurred, contrast-enhanced copy of C1
//       ("CLICK HERE") on the left, and the raw GFP channel on the right.
//    3. STAGE 1 -- CLICK: the multi-point tool is active. Click once in the
//       middle of every round / compact cell, then press OK. Each click is
//       grown into a cell outline with the wand (flood fill) tool.
//    4. STAGE 2 -- DRAW: the oval tool becomes active. For irregular cells,
//       hyphae, clumps, or anything the wand got wrong, draw an oval (or use
//       the freehand tool) around it and press "t" to add it. Press OK when done.
//    5. Review the ROI list, then Accept / Redo (with a new tolerance) / Skip / Quit.
//
//  OUTPUTS (per accepted image, written to the output folder you choose):
//    <name>_ROIs.zip   -- the ROI set, reopenable in Fiji at any time
//    <name>_ROIs.png   -- QC snapshot of C1 with the outlines drawn on it
//    interactive_cell_measurements.csv -- ONE master table, one row PER CELL,
//                                         appended as you go (crash-safe)
//
//  NOTE: because every cell is its own ROI, you get per-cell GFP values rather
//  than a single %Area per image, so no separate merge step is needed.
// ============================================================================

// ---------------- Defaults (all of these are editable in the startup dialog) --
var fileExtension        = ".czi";
var processSubfolders    = false; // true = walk into strain/condition subfolders too
var skipIfDone           = true;  // skip files that already have a saved _ROIs.zip

// Resume support. The file list is built in the same order every run (alphabetical
// within each folder), so the [n/total] numbers printed in the Log are stable as
// long as the input folder, extension and condition filter are unchanged.
var startAtImage         = 1;     // 1 = from the beginning; 104 = resume at image 104
var listFilesAndStop     = false; // print the numbered file list and quit, to check numbering

// Only process images whose CONDITION FOLDER name starts with this prefix.
// With "lc_" that means lc_t0, lc_t3, lc_t5, lc_t7, lc_t9.5, lc_t24 are processed
// and agar / glass / glass_pbs / plug / tapelift are ignored. Case-insensitive.
// Set useConditionFilter to false to process every condition.
var conditionFilter      = "lc_";
var useConditionFilter   = true;

var gfpChannelPrefix     = "C2";  // GFP channel (confirmed for this dataset)
var brightfieldPrefix    = "C1";  // brightfield channel used for selecting cells

var clickBlurSigma       = 2;     // blur applied ONLY to the click/wand copy of C1
var useAutoTolerance     = true;  // wand tolerance = factor x image SD (bit-depth independent)
var autoToleranceFactor  = 1.5;   // starting factor, and the fallback if you skip calibration
var manualTolerance      = 1000;  // used when useAutoTolerance is false

// ---- Tolerance calibration -------------------------------------------------
// The first image of each new condition folder shows the same colour-coded
// tolerance sweep as the preview macro: you click one typical cell, compare the
// outlines, and pick a factor. That factor is then reused for every remaining
// image in the folder -- multiplied by EACH image's own SD, so it still adapts
// to frame-by-frame contrast changes within the condition.
//   "Each folder"          -- recalibrate per strain+condition folder (e.g. 2.1/lc_t3)
//   "Each condition name"  -- one calibration per condition, shared across strains
//   "Once at start"        -- calibrate on the very first image only
//   "Never"                -- no calibration, just use the auto factor above
var calibrationScope     = "Each folder";
var calibrationFactors   = newArray(0.5, 1, 1.5, 2, 3);
var calibrationColors    = newArray("red", "orange", "yellow", "green", "cyan");
var currentFactor        = 1.5;   // set from autoToleranceFactor after the dialog
var doneKeys             = "|";   // calibration keys already handled this run
var wandMode             = "4-connected smooth";
var fillRoiHoles         = true;  // fill interiors so cells are solid, not rings
var smoothRoi            = true;  // interpolate the wand outline to remove jagged pixels

var minCellArea          = 20;    // reject wand results smaller than this
var maxCellAreaFraction  = 0.25;  // reject wand results bigger than 25% of the frame (runaway flood fill)

var measureInPixels      = true;  // ignore .czi micron calibration, report pixels
var gfpThreshold         = 0;     // >0 also reports % of pixels above this value per cell
var saveQcSnapshot       = true;

// ---------------- Startup dialog ------------------------------------------
Dialog.create("Interactive cell segmentation");
Dialog.addString("File extension:", fileExtension);
Dialog.addCheckbox("Include subfolders", processSubfolders);
Dialog.addCheckbox("Skip files that already have saved ROIs", skipIfDone);
Dialog.addNumber("Start at image number:", startAtImage);
Dialog.addCheckbox("Just list the numbered file order and stop", listFilesAndStop);
Dialog.addCheckbox("Only process conditions starting with:", useConditionFilter);
Dialog.addString("   condition prefix:", conditionFilter);
Dialog.addMessage("--- Click-to-select (wand) settings ---");
Dialog.addNumber("Blur sigma on click copy:", clickBlurSigma);
Dialog.addCheckbox("Auto wand tolerance (factor x image SD)", useAutoTolerance);
Dialog.addNumber("   starting / fallback factor:", autoToleranceFactor);
Dialog.addNumber("   manual tolerance (if auto off):", manualTolerance);
Dialog.addChoice("Calibrate tolerance:", newArray("Each folder", "Each condition name", "Once at start", "Never"), calibrationScope);
Dialog.addCheckbox("Fill holes inside each cell", fillRoiHoles);
Dialog.addCheckbox("Smooth cell outlines", smoothRoi);
Dialog.addNumber("Min cell area (px^2):", minCellArea);
Dialog.addMessage("--- Measurement ---");
Dialog.addCheckbox("Measure in pixels (ignore micron calibration)", measureInPixels);
Dialog.addNumber("GFP threshold (0 = off):", gfpThreshold);
Dialog.addCheckbox("Save QC snapshot PNG", saveQcSnapshot);
Dialog.show();

fileExtension       = Dialog.getString();
processSubfolders   = Dialog.getCheckbox();
skipIfDone          = Dialog.getCheckbox();
startAtImage        = Dialog.getNumber();
listFilesAndStop    = Dialog.getCheckbox();
useConditionFilter  = Dialog.getCheckbox();
conditionFilter     = Dialog.getString();
clickBlurSigma      = Dialog.getNumber();
useAutoTolerance    = Dialog.getCheckbox();
autoToleranceFactor = Dialog.getNumber();
manualTolerance     = Dialog.getNumber();
calibrationScope    = Dialog.getChoice();
currentFactor       = autoToleranceFactor;
fillRoiHoles        = Dialog.getCheckbox();
smoothRoi           = Dialog.getCheckbox();
minCellArea         = Dialog.getNumber();
measureInPixels     = Dialog.getCheckbox();
gfpThreshold        = Dialog.getNumber();
saveQcSnapshot      = Dialog.getCheckbox();

inputDir  = getDirectory("Choose the folder of .czi images to segment");
outputDir = getDirectory("Choose an output folder for ROIs and measurements");

csvPath = outputDir + "interactive_cell_measurements.csv";
if (!File.exists(csvPath)) {
    File.append("File,Folder,Strain,Condition,CellID,SelectionMethod,Area,GFP_Mean,GFP_Median,GFP_Max,GFP_RawIntDen,GFP_PctAreaAboveThreshold,GFP_Threshold,Units", csvPath);
}

// Record of which tolerance factor was chosen where, so the run is reproducible
calibPath = outputDir + "calibration_log.csv";
if (!File.exists(calibPath)) {
    File.append("CalibrationKey,CalibratedOnFile,Strain,Condition,ChosenFactor,ImageSD,ResultingTolerance", calibPath);
}

// ---------------- Helpers --------------------------------------------------
var collected = newArray(0);
var excludedCount = 0;

// A file passes if the folder it sits in (its condition folder) starts with the prefix
function conditionAllowed(path) {
    if (!useConditionFilter) return true;
    if (conditionFilter == "") return true;
    cond = File.getName(File.getParent(path));
    return startsWith(toLowerCase(cond), toLowerCase(conditionFilter));
}

function collectFiles(dir) {
    list = getFileList(dir);
    for (i = 0; i < list.length; i++) {
        if (endsWith(list[i], "/")) {
            if (processSubfolders) collectFiles(dir + list[i]);
        } else if (endsWith(toLowerCase(list[i]), toLowerCase(fileExtension))) {
            if (conditionAllowed(dir + list[i])) {
                collected = Array.concat(collected, dir + list[i]);
            } else {
                excludedCount++;
            }
        }
    }
}

// If a channel window is a Z-stack, keep only the middle slice
// ---- Opening .czi files ----------------------------------------------------
// "Bio-Formats Windowless Importer" applies the options string WITHOUT showing
// the import dialog. The plain "Bio-Formats Importer" can pop that dialog on a
// thread the macro cannot answer, which looks exactly like the macro hanging or
// timing out. Change this back only if you need the dialog.
var importerCommand = "Bio-Formats Windowless Importer";

var openedC1 = "";
var openedC2 = "";
var openedW  = 0;
var openedH  = 0;

// Opens a file and leaves the brightfield and GFP channels in openedC1/openedC2.
// Copes with the windowless importer having "split channels" saved in its own
// preferences, in which case the channels arrive already separated.
function openAndSplitChannels(path) {
    openedC1 = "";
    openedC2 = "";
    openedW  = 0;
    openedH  = 0;

    before = getList("image.titles");
    run(importerCommand, "open=[" + path + "] autoscale color_mode=Default view=Hyperstack stack_order=XYCZT");
    after = getList("image.titles");

    fresh = newArray(0);
    for (i = 0; i < after.length; i++) {
        wasOpen = false;
        for (j = 0; j < before.length; j++) { if (after[i] == before[j]) wasOpen = true; }
        if (!wasOpen) fresh = Array.concat(fresh, after[i]);
    }
    if (fresh.length == 0) {
        print("   nothing opened -- Bio-Formats could not read this file.");
        return;
    }

    selectWindow(fresh[0]);
    getDimensions(w0, h0, ch, sl, fr);
    openedW = w0;
    openedH = h0;

    // One window with several channels -> split it ourselves
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
        print("   only one channel in this file -- skipping.");
        selectWindow(fresh[0]);
        close();
        return;
    }

    // Identify channels by the C1-/C2- prefix Bio-Formats and Split Channels use
    leftovers = newArray(0);
    for (i = 0; i < fresh.length; i++) {
        if (startsWith(fresh[i], brightfieldPrefix + "-") || startsWith(fresh[i], brightfieldPrefix + "_")) {
            openedC1 = fresh[i];
        } else if (startsWith(fresh[i], gfpChannelPrefix + "-") || startsWith(fresh[i], gfpChannelPrefix + "_")) {
            openedC2 = fresh[i];
        } else {
            leftovers = Array.concat(leftovers, fresh[i]);
        }
    }

    // Fallback: no C1-/C2- prefixes, but exactly two windows appeared. Assume
    // channel order (first = brightfield, second = GFP) and say so in the Log.
    if (openedC1 == "" && openedC2 == "" && leftovers.length == 2) {
        openedC1 = leftovers[0];
        openedC2 = leftovers[1];
        leftovers = newArray(0);
        print("   note: no C1-/C2- prefixes found; using window order "
              + "(1st = brightfield, 2nd = GFP). Check this is right.");
    }

    for (i = 0; i < leftovers.length; i++) {
        if (isOpen(leftovers[i])) { selectWindow(leftovers[i]); close(); }
    }
}

// Which images share one calibration, per the chosen scope
function calibrationKey(path) {
    if (calibrationScope == "Never")             return "";
    if (calibrationScope == "Once at start")     return "ALL";
    if (calibrationScope == "Each condition name") return File.getName(File.getParent(path));
    return File.getParent(path); // "Each folder" -- full strain+condition path
}

// Show what the wand grabs at several tolerances at once, and let the user pick.
// Returns the chosen factor (a multiple of this image's SD).
function calibrateTolerance(imgTitle, imgSD, maxArea, headerText) {
    factor = currentFactor;
    calibrating = true;
    while (calibrating) {
        calibrating = false;
        selectWindow(imgTitle);
        Overlay.remove;
        run("Select None");
        setTool("point");
        waitForUser("TOLERANCE CALIBRATION\n \n" + headerText + "\n \n"
            + "Click ONCE in the middle of a typical, well-separated cell,\n"
            + "then press OK.\n \n"
            + "The macro will draw what the wand grabs at five tolerances so\n"
            + "you can pick the one that traces the real cell edge. This\n"
            + "choice is then reused for the rest of this condition.");

        selectWindow(imgTitle);
        cx = newArray(0);
        cy = newArray(0);
        if (selectionType() == 10 || selectionType() == 0 || selectionType() == 7) {
            getSelectionCoordinates(cx, cy);
        }
        run("Select None");

        if (cx.length == 0) {
            print("   no calibration point clicked -- keeping factor " + factor);
            Overlay.remove;
            return factor;
        }

        px = cx[cx.length - 1];
        py = cy[cy.length - 1];
        labels = newArray(calibrationFactors.length);
        print("   calibration sweep at (" + px + "," + py + "), image SD = " + d2s(imgSD, 1) + ":");
        for (t = 0; t < calibrationFactors.length; t++) {
            tol = calibrationFactors[t] * imgSD;
            if (tol < 1) tol = 1;
            selectWindow(imgTitle);
            run("Select None");
            doWand(px, py, tol, wandMode);
            if (selectionType() == -1) {
                labels[t] = "" + calibrationFactors[t] + "  (" + calibrationColors[t] + ")  nothing found";
                print("      " + calibrationFactors[t] + " x SD = " + d2s(tol, 1) + " -> nothing");
            } else {
                aa = getValue("Area");
                note = "";
                if (aa > maxArea) note = "  TOO BIG - flooded";
                else if (aa < minCellArea) note = "  too small";
                Overlay.addSelection(calibrationColors[t], 2);
                labels[t] = "" + calibrationFactors[t] + "  (" + calibrationColors[t] + ")  area " + d2s(aa, 0) + " px" + note;
                print("      " + calibrationFactors[t] + " x SD = " + d2s(tol, 1)
                      + " -> area " + d2s(aa, 0) + " px  (" + calibrationColors[t] + ")" + note);
            }
        }
        selectWindow(imgTitle);
        run("Select None");
        Overlay.show;

        Dialog.create("Pick a tolerance for this condition");
        Dialog.addMessage(headerText + "\n \nCompare the coloured outlines.\n"
            + "Prefer a factor where the area barely changes from its neighbours\n"
            + "-- that plateau means the wand is stopping at a real edge.");
        Dialog.addRadioButtonGroup("Tolerance (x image SD):", labels, labels.length, 1, labels[2]);
        Dialog.addCheckbox("Calibrate on a different cell instead", false);
        Dialog.show();
        pick  = Dialog.getRadioButton();
        again = Dialog.getCheckbox();

        for (t = 0; t < labels.length; t++) {
            if (labels[t] == pick) factor = calibrationFactors[t];
        }
        if (again) calibrating = true;
    }
    selectWindow(imgTitle);
    Overlay.remove;
    run("Select None");
    return factor;
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

// Turn the current selection into a hole-free version of itself
// Fill interior holes so a cell is solid, not a ring.
// Create Mask / Create Selection can invert depending on the binary LUT and the
// black-background setting, which would return the BACKGROUND instead of the
// cell. Normalising the LUT and selecting via an explicit threshold avoids that.
function fillHolesInSelection(hostTitle) {
    if (selectionType() == -1) return;
    selectWindow(hostTitle);
    run("Create Mask");
    maskTitle = getTitle();
    if (is("Inverting LUT")) run("Invert LUT");   // cell = white(255), background = black(0)
    setOption("BlackBackground", true);
    run("Fill Holes");
    setThreshold(128, 255);                        // select white pixels explicitly
    run("Create Selection");
    resetThreshold();
    if (selectionType() == -1) {
        close(maskTitle);
        selectWindow(hostTitle);
        return;
    }
    selectWindow(hostTitle);
    run("Restore Selection");
    close(maskTitle);
    selectWindow(hostTitle);
}

function upperLimitForBitDepth() {
    bd = bitDepth();
    if (bd == 8)  return 255;
    if (bd == 16) return 65535;
    return 1e30; // 32-bit float
}

function safeClose(t) {
    if (isOpen(t)) { selectWindow(t); close(); }
}

// ---------------- Setup ----------------------------------------------------
setBatchMode(false);           // interaction requires visible windows
setOption("BlackBackground", true);
run("Set Measurements...", "area mean min median integrated area_fraction display redirect=None decimal=4");
roiManager("reset");
run("Clear Results");
if (measureInPixels) units = "pixels"; else units = "calibrated";

collectFiles(inputDir);
print("\\Clear");
if (collected.length == 0) {
    msg = "No " + fileExtension + " files to process in " + inputDir;
    if (useConditionFilter) {
        msg = msg + "\n \nThe condition filter \"" + conditionFilter + "\" excluded "
            + excludedCount + " file(s).\nIf your images sit in strain folders "
            + "(1.1/2.1/3.1), tick 'Include subfolders'\nand choose the folder above them.";
    }
    exit(msg);
}
print("Found " + collected.length + " " + fileExtension + " file(s) to process.");
if (useConditionFilter) {
    print("Condition filter \"" + conditionFilter + "\" excluded " + excludedCount + " file(s) in other conditions.");
}

if (listFilesAndStop) {
    print("");
    print("--- File order (use these numbers for 'Start at image number') ---");
    for (i = 0; i < collected.length; i++) {
        marker = "";
        rz = outputDir + replace(File.getName(collected[i]), fileExtension, "") + "_ROIs.zip";
        if (File.exists(rz)) marker = "   [already has ROIs]";
        print("  " + (i + 1) + ": " + File.getName(File.getParent(collected[i]))
              + "/" + File.getName(collected[i]) + marker);
    }
    exit("Listed " + collected.length + " file(s) in the Log. Nothing was processed.");
}

if (startAtImage > 1) {
    if (startAtImage > collected.length)
        exit("Start at image number (" + startAtImage + ") is beyond the "
             + collected.length + " matched file(s).");
    print("Resuming at image " + startAtImage + ": "
          + File.getName(collected[startAtImage - 1])
          + "  (skipping the first " + (startAtImage - 1) + ")");
}

quit = false;

// ---------------- Main loop ------------------------------------------------
for (f = 0; f < collected.length && !quit; f++) {

    if (f + 1 < startAtImage) continue;   // resume point

    filePath  = collected[f];
    fileName  = File.getName(filePath);
    shortName = replace(fileName, fileExtension, "");
    folder    = File.getParent(filePath);
    condition = File.getName(folder);
    strain    = File.getName(File.getParent(folder));

    roiZipPath = outputDir + shortName + "_ROIs.zip";
    if (skipIfDone && File.exists(roiZipPath)) {
        print("[" + (f + 1) + "/" + collected.length + "] " + fileName + " -- already has ROIs, skipping.");
        continue;
    }

    imgTolOverride = 0;   // set only when you choose "Redo this image" with a typed tolerance
    forceRecalibrate = false;

    redo = true;
    while (redo && !quit) {
        redo = false;
        roiManager("reset");

        // ---- Open and split ---------------------------------------------
        openAndSplitChannels(filePath);
        c1Title = openedC1;
        c2Title = openedC2;
        w = openedW;
        h = openedH;
        if (c1Title == "" || c2Title == "") {
            print("Skipping " + fileName + " -- could not get both channels.");
            close("*");
            continue;
        }

        reduceToMiddleSlice(c1Title);
        reduceToMiddleSlice(c2Title);

        if (measureInPixels) {
            selectWindow(c1Title); run("Set Scale...", "distance=0 known=0 unit=pixel");
            selectWindow(c2Title); run("Set Scale...", "distance=0 known=0 unit=pixel");
        }

        // ---- Build the interactive click target ----------------------------
        // A blurred, contrast-enhanced copy of C1. The blur stops the wand from
        // stopping at internal organelle texture; Enhance Contrast is display-only.
        clickTitle = shortName + " -- CLICK HERE (C1)";
        selectWindow(c1Title);
        run("Duplicate...", "title=[" + clickTitle + "] duplicate");
        if (clickBlurSigma > 0) run("Gaussian Blur...", "sigma=" + clickBlurSigma);
        run("Enhance Contrast", "saturated=0.35");   // display only, pixels unchanged

        gfpViewTitle = shortName + " -- GFP (C2)";
        selectWindow(c2Title);
        run("Duplicate...", "title=[" + gfpViewTitle + "] duplicate");
        run("Enhance Contrast", "saturated=0.35");   // display only
        setLocation(560, 100);

        selectWindow(clickTitle);
        setLocation(40, 100);

        // ---- Wand tolerance -------------------------------------------------
        getStatistics(statArea, statMean, statMin, statMax, statSD);
        maxCellArea = maxCellAreaFraction * w * h;

        // Calibrate on the first image of each new condition folder (per the scope
        // setting), then reuse that factor for the rest of the folder.
        calKey = calibrationKey(filePath);
        if (calKey != "" && useAutoTolerance) {
            if (forceRecalibrate || indexOf(doneKeys, "|" + calKey + "|") < 0) {
                print("--- Calibrating for: " + strain + " / " + condition + " ---");
                currentFactor = calibrateTolerance(clickTitle, statSD, maxCellArea,
                    strain + " / " + condition + "\n" + fileName);
                if (indexOf(doneKeys, "|" + calKey + "|") < 0) doneKeys = doneKeys + calKey + "|";
                forceRecalibrate = false;
                File.append(calKey + "," + fileName + "," + strain + "," + condition + ","
                    + currentFactor + "," + d2s(statSD, 2) + "," + d2s(currentFactor * statSD, 2), calibPath);
                print("   chosen factor for this condition: " + currentFactor + " x SD");
            }
        }

        if (imgTolOverride > 0) {
            tolerance = imgTolOverride;   // this image only; does not change the folder's factor
        } else if (useAutoTolerance) {
            tolerance = currentFactor * statSD;
        } else {
            tolerance = manualTolerance;
        }
        if (tolerance < 1) tolerance = 1;

        roiManager("reset");
        roiManager("Show All with labels");

        // ---- STAGE 1: click each round cell --------------------------------
        setTool("multipoint");
        selectWindow(clickTitle);
        run("Select None");
        waitForUser("STAGE 1 of 2  --  CLICK CELLS\n \n"
            + fileName + "   (" + (f + 1) + " of " + collected.length + ")\n \n"
            + "Click once in the MIDDLE of every round / compact cell in the\n"
            + "'CLICK HERE' window. Shift-click is not needed - the multi-point\n"
            + "tool keeps every click.\n \n"
            + "Mis-clicked? Alt-click the point to remove it.\n"
            + "No round cells in this image? Just press OK.\n \n"
            + "Wand tolerance for this image: " + d2s(tolerance, 1));

        clickX = newArray(0);
        clickY = newArray(0);
        if (isOpen(clickTitle)) {
            selectWindow(clickTitle);
            if (selectionType() == 10 || selectionType() == 0 || selectionType() == 7) {
                getSelectionCoordinates(clickX, clickY);
            }
        }
        run("Select None");

        accepted = 0;
        rejected = 0;
        for (p = 0; p < clickX.length; p++) {
            selectWindow(clickTitle);
            run("Select None");
            doWand(clickX[p], clickY[p], tolerance, wandMode);
            if (selectionType() == -1) {
                rejected++;
                print("   click " + (p + 1) + " rejected: the wand found nothing at that point");
                continue;
            }
            aRaw = getValue("Area");

            if (fillRoiHoles) {
                fillHolesInSelection(clickTitle);
                aFilled = 0;
                if (selectionType() != -1) aFilled = getValue("Area");
                // Safety net: if filling holes ballooned the area it has grabbed the
                // background, so fall back to the plain wand outline for this cell.
                if (aFilled <= 0 || aFilled > 3 * aRaw) {
                    selectWindow(clickTitle);
                    run("Select None");
                    doWand(clickX[p], clickY[p], tolerance, wandMode);
                    print("   click " + (p + 1) + ": fill-holes misfired (" + d2s(aRaw, 0)
                          + " -> " + d2s(aFilled, 0) + " px), kept the plain outline");
                }
            }
            if (smoothRoi && selectionType() != -1) run("Interpolate", "interval=2 smooth adjust");
            if (selectionType() == -1) { rejected++; continue; }
            a = getValue("Area");
            if (a < minCellArea) {
                rejected++;
                print("   click " + (p + 1) + " rejected: area " + d2s(a, 1)
                      + " px is under the " + minCellArea + " px minimum -- raise the tolerance");
                continue;
            }
            if (a > maxCellArea) {
                rejected++;
                print("   click " + (p + 1) + " rejected: area " + d2s(a, 1)
                      + " px is over the " + d2s(maxCellArea, 0)
                      + " px flood limit -- lower the tolerance");
                continue;
            }
            roiManager("Add");
            roiManager("Select", roiManager("count") - 1);
            roiManager("Rename", "click_" + IJ.pad(accepted + 1, 3));
            accepted++;
        }
        selectWindow(clickTitle);
        run("Select None");
        clickCount = roiManager("count");
        roiManager("Deselect");
        roiManager("Show All with labels");
        print("[" + (f + 1) + "/" + collected.length + "] " + fileName + ": "
              + accepted + " cell(s) from clicks, " + rejected + " click(s) rejected.");
        if (accepted == 0 && clickX.length > 0)
            print("   >> every click was rejected. See the reason above: 'under the minimum'"
                  + " means raise the tolerance, 'over the flood limit' means lower it.");

        // ---- STAGE 2: draw around irregular cells ---------------------------
        setTool("oval");
        selectWindow(clickTitle);
        waitForUser("STAGE 2 of 2  --  DRAW IRREGULAR CELLS\n \n"
            + "For any cell the wand missed or got wrong (elongated cells,\n"
            + "hyphae, clumps, touching cells):\n \n"
            + "   1. Drag an OVAL around it  (or pick the freehand tool for\n"
            + "      awkward shapes - both work).\n"
            + "   2. Press the 't' key to add it to the ROI Manager.\n"
            + "   3. Repeat for each one.\n \n"
            + "Delete a bad ROI by selecting it in the ROI Manager and\n"
            + "pressing Delete. Nothing to draw? Just press OK.\n \n"
            + "Currently in the list: " + clickCount + " cell(s).");

        selectWindow(clickTitle);
        run("Select None");
        totalRois = roiManager("count");
        for (r = clickCount; r < totalRois; r++) {
            roiManager("Select", r);
            roiManager("Rename", "manual_" + IJ.pad(r - clickCount + 1, 3));
        }
        roiManager("Deselect");

        // ---- Review ---------------------------------------------------------
        Dialog.create("Review  --  " + fileName);
        Dialog.addMessage("ROIs currently selected: " + roiManager("count")
            + "\n(" + clickCount + " clicked, " + (roiManager("count") - clickCount) + " drawn)");
        Dialog.addRadioButtonGroup("What next?",
            newArray("Accept and continue", "Redo this image", "Recalibrate this condition and redo",
                     "Skip this image (save nothing)", "Quit"),
            5, 1, "Accept and continue");
        Dialog.addNumber("Wand tolerance for 'Redo this image':", tolerance);
        Dialog.show();
        action        = Dialog.getRadioButton();
        redoTolerance = Dialog.getNumber();

        if (action == "Redo this image") {
            imgTolOverride = redoTolerance;   // applies to THIS image only
            redo = true;
        } else if (action == "Recalibrate this condition and redo") {
            imgTolOverride   = 0;             // go back to the factor-based tolerance
            forceRecalibrate = true;          // re-run the sweep, update the folder's factor
            redo = true;
        } else if (action == "Quit") {
            quit = true;
        } else if (action == "Accept and continue" && roiManager("count") > 0) {

            // ---- Save ROI set -----------------------------------------------
            roiManager("Deselect");
            roiManager("Save", roiZipPath);

            // ---- Measure every ROI on the RAW GFP channel --------------------
            upper = 0;
            selectWindow(c2Title);
            upper = upperLimitForBitDepth();
            n = roiManager("count");
            for (r = 0; r < n; r++) {
                selectWindow(c2Title);
                roiManager("Select", r);
                roiName = Roi.getName;
                if (startsWith(roiName, "manual")) method = "drawn"; else method = "clicked";

                resetThreshold();
                a      = getValue("Area");
                mean   = getValue("Mean");
                median = getValue("Median");
                mx     = getValue("Max");
                intden = getValue("RawIntDen");

                pctAbove = -1;
                if (gfpThreshold > 0) {
                    setThreshold(gfpThreshold, upper);
                    pctAbove = getValue("%Area");
                    resetThreshold();
                }

                row = fileName + "," + condition + "," + strain + "," + condition + ","
                    + roiName + "," + method + ","
                    + d2s(a, 4) + "," + d2s(mean, 4) + "," + d2s(median, 4) + ","
                    + d2s(mx, 4) + "," + d2s(intden, 4) + ","
                    + d2s(pctAbove, 4) + "," + gfpThreshold + "," + units;
                File.append(row, csvPath);
            }
            selectWindow(c2Title);
            run("Select None");
            resetThreshold();

            // ---- QC snapshot -------------------------------------------------
            if (saveQcSnapshot) {
                selectWindow(c1Title);
                run("Enhance Contrast", "saturated=0.35");
                roiManager("Deselect");
                roiManager("Show All without labels");
                run("From ROI Manager");
                run("Flatten");
                saveAs("PNG", outputDir + shortName + "_ROIs.png");
                close();
                selectWindow(c1Title);
                run("Remove Overlay");
            }
            print("   saved " + roiManager("count") + " ROI(s) -> " + shortName + "_ROIs.zip");
        } else if (action == "Accept and continue") {
            print("   no ROIs selected -- nothing saved for " + fileName);
        } else {
            print("   skipped " + fileName);
        }

        // ---- Clean up this image -------------------------------------------
        safeClose(clickTitle);
        safeClose(gfpViewTitle);
        safeClose(c1Title);
        safeClose(c2Title);
    }
}

roiManager("reset");
setTool("rectangle");
print("----------------------------------------------------------");
print("Done. Per-cell measurements: " + csvPath);
print("Tolerance factors used per condition: " + calibPath);
print("ROI sets (.zip) and QC snapshots (.png) are in: " + outputDir);
print("To revisit an image: open it, then ROI Manager > More > Open, and load its _ROIs.zip.");
