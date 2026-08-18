// ============================================================================
//  PREVIEW: Interactive Cell Segmentation  --  GFP quantification project
// ============================================================================
//  Test-drive version of interactive_cell_segmentation.ijm. Same click-to-select
//  workflow, but on only a few images and it SAVES NOTHING -- no ROI zips, no
//  CSV, no PNGs. Use it to find a wand tolerance that works for your cell types,
//  then put that number into the full macro.
//
//  WHAT IT ADDS OVER THE FULL MACRO:
//    STAGE 0 -- TOLERANCE CALIBRATION. You click ONE typical cell and the macro
//    draws what the wand would grab at five different tolerances, colour-coded,
//    all at once. You pick the one that traces the cell edge properly. No more
//    guessing at a number that means nothing until you see it.
//
//  Measurements go to the Results window on screen so you can eyeball whether
//  the per-cell GFP values look sensible. Windows are left open at the end with
//  the outlines overlaid on both the raw brightfield and the raw GFP channel.
// ============================================================================

// ---------------- Defaults (editable in the startup dialog) -----------------
var fileExtension        = ".czi";
var processSubfolders    = false;
var maxPreviewImages     = 2;     // how many files to try

// Only preview images whose CONDITION FOLDER starts with this (case-insensitive).
// "lc_" = lc_t0, lc_t3, lc_t5, lc_t7, lc_t9.5, lc_t24 only.
var conditionFilter      = "lc_";
var useConditionFilter   = true;

var gfpChannelPrefix     = "C2";
var brightfieldPrefix    = "C1";

var clickBlurSigma       = 2;
var wandMode             = "4-connected smooth";
var fillRoiHoles         = true;
var smoothRoi            = true;
var minCellArea          = 20;
var maxCellAreaFraction  = 0.25;

var measureInPixels      = true;
var gfpThreshold         = 0;     // >0 also reports % of pixels above this value per cell

// Tolerances offered in the calibration step, as multiples of the image SD
var calibrationFactors   = newArray(0.5, 1, 1.5, 2, 3);
var calibrationColors    = newArray("red", "orange", "yellow", "green", "cyan");

// ---------------- Startup dialog -------------------------------------------
Dialog.create("PREVIEW  --  interactive cell segmentation");
Dialog.addMessage("Nothing is saved. This is for tuning only.");
Dialog.addString("File extension:", fileExtension);
Dialog.addNumber("How many images to preview:", maxPreviewImages);
Dialog.addCheckbox("Include subfolders", processSubfolders);
Dialog.addCheckbox("Only preview conditions starting with:", useConditionFilter);
Dialog.addString("   condition prefix:", conditionFilter);
Dialog.addNumber("Blur sigma on click copy:", clickBlurSigma);
Dialog.addCheckbox("Fill holes inside each cell", fillRoiHoles);
Dialog.addCheckbox("Smooth cell outlines", smoothRoi);
Dialog.addNumber("Min cell area (px^2):", minCellArea);
Dialog.addCheckbox("Measure in pixels (ignore micron calibration)", measureInPixels);
Dialog.addNumber("GFP threshold (0 = off):", gfpThreshold);
Dialog.show();

fileExtension      = Dialog.getString();
maxPreviewImages   = Dialog.getNumber();
processSubfolders  = Dialog.getCheckbox();
useConditionFilter = Dialog.getCheckbox();
conditionFilter    = Dialog.getString();
clickBlurSigma     = Dialog.getNumber();
fillRoiHoles       = Dialog.getCheckbox();
smoothRoi          = Dialog.getCheckbox();
minCellArea        = Dialog.getNumber();
measureInPixels    = Dialog.getCheckbox();
gfpThreshold       = Dialog.getNumber();

inputDir = getDirectory("Choose a folder of test images (an lc_ condition folder, or its parent with subfolders on)");

// ---------------- Helpers ---------------------------------------------------
var collected = newArray(0);
var excludedCount = 0;

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
            if (conditionAllowed(dir + list[i])) collected = Array.concat(collected, dir + list[i]);
            else excludedCount++;
        }
    }
}

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
    return 1e30;
}

function safeClose(t) {
    if (isOpen(t)) { selectWindow(t); close(); }
}

// ---------------- Setup ------------------------------------------------------
setBatchMode(false);
setOption("BlackBackground", true);
run("Set Measurements...", "area mean min median integrated area_fraction display redirect=None decimal=4");
roiManager("reset");
run("Clear Results");
if (measureInPixels) units = "pixels"; else units = "calibrated";

collectFiles(inputDir);
print("\\Clear");
if (collected.length == 0) {
    msg = "No " + fileExtension + " files to preview in " + inputDir;
    if (useConditionFilter) {
        msg = msg + "\n \nThe condition filter \"" + conditionFilter + "\" excluded "
            + excludedCount + " file(s).\nIf your images sit in strain folders "
            + "(1.1/2.1/3.1), tick 'Include subfolders'\nand choose the folder above them.";
    }
    exit(msg);
}
print("PREVIEW MODE -- nothing will be saved.");
print("Matched " + collected.length + " file(s); previewing up to " + maxPreviewImages + ".");
if (useConditionFilter) print("Condition filter \"" + conditionFilter + "\" excluded " + excludedCount + " file(s).");
print("");

previewed = 0;
quit = false;

// ---------------- Main loop --------------------------------------------------
for (f = 0; f < collected.length && previewed < maxPreviewImages && !quit; f++) {

    filePath  = collected[f];
    fileName  = File.getName(filePath);
    shortName = replace(fileName, fileExtension, "");
    condition = File.getName(File.getParent(filePath));
    strain    = File.getName(File.getParent(File.getParent(filePath)));

    retry = true;
    while (retry && !quit) {
        retry = false;
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

        // ---- Windows -----------------------------------------------------
        clickTitle = shortName + " -- CLICK HERE (C1)";
        selectWindow(c1Title);
        run("Duplicate...", "title=[" + clickTitle + "] duplicate");
        if (clickBlurSigma > 0) run("Gaussian Blur...", "sigma=" + clickBlurSigma);
        run("Enhance Contrast", "saturated=0.35");   // display only

        gfpViewTitle = shortName + " -- GFP (C2)";
        selectWindow(c2Title);
        run("Duplicate...", "title=[" + gfpViewTitle + "] duplicate");
        run("Enhance Contrast", "saturated=0.35");
        setLocation(560, 100);

        selectWindow(clickTitle);
        setLocation(40, 100);
        getStatistics(statArea, statMean, statMin, statMax, statSD);
        maxCellArea = maxCellAreaFraction * w * h;

        print("=== " + strain + " / " + condition + " / " + fileName + " ===");
        print("   image SD = " + d2s(statSD, 1) + ", bit depth = " + bitDepth());

        // ---- STAGE 0: tolerance calibration --------------------------------
        tolerance = 1.5 * statSD;
        chosenFactor = 1.5;
        calibrating = true;
        while (calibrating) {
            calibrating = false;
            selectWindow(clickTitle);
            Overlay.remove;
            run("Select None");
            setTool("point");
            waitForUser("STAGE 0  --  TOLERANCE CALIBRATION\n \n"
                + fileName + "\n \n"
                + "Click ONCE in the middle of a typical, well-separated cell,\n"
                + "then press OK.\n \n"
                + "The macro will draw what the wand grabs at five different\n"
                + "tolerances so you can pick the one that traces the real edge.");

            selectWindow(clickTitle);
            calX = newArray(0);
            calY = newArray(0);
            if (selectionType() == 10 || selectionType() == 0 || selectionType() == 7) {
                getSelectionCoordinates(calX, calY);
            }
            run("Select None");
            if (calX.length == 0) {
                print("   no calibration point clicked -- using 1.5 x SD = " + d2s(tolerance, 1));
            } else {
                px = calX[calX.length - 1];
                py = calY[calY.length - 1];
                labels = newArray(calibrationFactors.length);
                print("   calibration at (" + px + "," + py + "):");
                for (t = 0; t < calibrationFactors.length; t++) {
                    tol = calibrationFactors[t] * statSD;
                    if (tol < 1) tol = 1;
                    selectWindow(clickTitle);
                    run("Select None");
                    doWand(px, py, tol, wandMode);
                    if (selectionType() == -1) {
                        labels[t] = "" + calibrationFactors[t] + "  (" + calibrationColors[t] + ")  nothing found";
                        print("      " + calibrationFactors[t] + " x SD = " + d2s(tol, 1) + " -> nothing (" + calibrationColors[t] + ")");
                    } else {
                        aa = getValue("Area");
                        note = "";
                        if (aa > maxCellArea) note = "  TOO BIG - flooded";
                        else if (aa < minCellArea) note = "  too small";
                        Overlay.addSelection(calibrationColors[t], 2);
                        labels[t] = "" + calibrationFactors[t] + "  (" + calibrationColors[t] + ")  area " + d2s(aa, 0) + " px" + note;
                        print("      " + calibrationFactors[t] + " x SD = " + d2s(tol, 1)
                              + " -> area " + d2s(aa, 0) + " px  (" + calibrationColors[t] + ")" + note);
                    }
                }
                selectWindow(clickTitle);
                run("Select None");
                Overlay.show;

                Dialog.create("Pick a tolerance");
                Dialog.addMessage("Compare the coloured outlines on the image.\nPick the one that follows the cell edge:");
                Dialog.addRadioButtonGroup("Tolerance (x image SD):", labels, labels.length, 1, labels[2]);
                Dialog.addCheckbox("Actually, let me calibrate on a different cell", false);
                Dialog.show();
                pick = Dialog.getRadioButton();
                again = Dialog.getCheckbox();

                for (t = 0; t < labels.length; t++) {
                    if (labels[t] == pick) chosenFactor = calibrationFactors[t];
                }
                tolerance = chosenFactor * statSD;
                if (tolerance < 1) tolerance = 1;
                if (again) calibrating = true;
            }
        }
        selectWindow(clickTitle);
        Overlay.remove;
        print("   using tolerance " + d2s(tolerance, 1) + "  (= " + chosenFactor + " x SD)");

        roiManager("reset");
        roiManager("Show All with labels");

        // ---- STAGE 1: click round cells -----------------------------------
        setTool("multipoint");
        selectWindow(clickTitle);
        run("Select None");
        waitForUser("STAGE 1 of 2  --  CLICK CELLS\n \n"
            + "Click once in the MIDDLE of every round / compact cell,\n"
            + "then press OK. Alt-click a point to remove it.\n \n"
            + "Tolerance in use: " + d2s(tolerance, 1));

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
        print("   " + accepted + " cell(s) from clicks, " + rejected + " rejected.");
        if (accepted == 0 && clickX.length > 0)
            print("   >> every click was rejected. See the reason above: 'under the minimum'"
                  + " means raise the tolerance, 'over the flood limit' means lower it.");

        // ---- STAGE 2: draw irregular cells ---------------------------------
        setTool("oval");
        selectWindow(clickTitle);
        waitForUser("STAGE 2 of 2  --  DRAW IRREGULAR CELLS\n \n"
            + "Drag an OVAL around anything the wand missed or got wrong\n"
            + "(elongated cells, hyphae, clumps), and press 't' to add it.\n"
            + "Freehand tool works too. Repeat for each one.\n \n"
            + "Currently in the list: " + clickCount + " cell(s).");

        selectWindow(clickTitle);
        run("Select None");
        totalRois = roiManager("count");
        for (r = clickCount; r < totalRois; r++) {
            roiManager("Select", r);
            roiManager("Rename", "manual_" + IJ.pad(r - clickCount + 1, 3));
        }
        roiManager("Deselect");

        // ---- Measure into the on-screen Results table (nothing saved) -------
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

            row = nResults;
            setResult("Label", row, fileName);
            setResult("Strain", row, strain);
            setResult("Condition", row, condition);
            setResult("CellID", row, roiName);
            setResult("Method", row, method);
            setResult("Area", row, a);
            setResult("GFP_Mean", row, mean);
            setResult("GFP_Median", row, median);
            setResult("GFP_Max", row, mx);
            setResult("GFP_RawIntDen", row, intden);
            if (gfpThreshold > 0) setResult("GFP_PctAbove", row, pctAbove);
            setResult("Units", row, units);
        }
        updateResults();
        selectWindow(c2Title);
        run("Select None");
        resetThreshold();

        // ---- Bake the outlines onto the raw channels for inspection ---------
        if (n > 0) {
            roiManager("Deselect");
            roiManager("Show All without labels");
            selectWindow(c1Title);
            rename(shortName + "_C1_raw+ROIs");
            c1Title = getTitle();
            run("Enhance Contrast", "saturated=0.35");
            run("From ROI Manager");
            setLocation(40, 100);

            selectWindow(c2Title);
            rename(shortName + "_C2_raw+ROIs");
            c2Title = getTitle();
            run("Enhance Contrast", "saturated=0.35");
            run("From ROI Manager");
            setLocation(560, 100);
        }
        safeClose(clickTitle);
        safeClose(gfpViewTitle);

        // ---- What next? -----------------------------------------------------
        Dialog.create("Preview  --  " + fileName);
        Dialog.addMessage("Selected " + n + " cell(s).\n"
            + "Check the outlines on the two raw windows and the Results table.\n"
            + "Nothing has been saved.");
        Dialog.addRadioButtonGroup("What next?",
            newArray("Next image", "Try this image again", "Stop here"), 3, 1, "Next image");
        Dialog.addCheckbox("Close this image's windows before continuing", true);
        Dialog.show();
        action    = Dialog.getRadioButton();
        closeThem = Dialog.getCheckbox();

        if (closeThem) {
            safeClose(c1Title);
            safeClose(c2Title);
        } else {
            selectWindow(c1Title); run("Select None");
            selectWindow(c2Title); run("Select None");
        }

        if (action == "Try this image again") {
            retry = true;
            if (closeThem == false) { safeClose(c1Title); safeClose(c2Title); }
        } else if (action == "Stop here") {
            quit = true;
        } else {
            previewed++;
        }
        print("");
    }
}

roiManager("reset");
setTool("rectangle");
print("----------------------------------------------------------");
print("Preview finished. Nothing was saved.");
print("Note the tolerance factor that worked (printed above for each image).");
print("In interactive_cell_segmentation.ijm, either leave 'Auto wand tolerance'");
print("ticked and set the auto factor to that number, or untick it and enter the");
print("raw tolerance value directly.");
