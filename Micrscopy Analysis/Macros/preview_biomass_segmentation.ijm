// ==== Preview: Biomass Mask Segmentation Check ====
// Run this BEFORE the full batch macro to confirm the brightfield (C1)
// segmentation correctly outlines fungal material (whether discrete cells
// or mycelium), not empty background.
//
// For a few .czi files, this opens the raw C1 (brightfield) and C2 (GFP)
// channels, computes the biomass mask from C1, and overlays its outline on
// BOTH channels so you can visually judge whether the mask is capturing the
// right thing. Nothing is measured, thresholded for GFP, or saved -- this is
// purely a visual sanity check. All windows are left open for you to inspect.
//
// If the mask looks INVERTED (outlining empty space instead of fungal
// material), change biomassDarkBackground below and re-run this preview.

fileExtension = ".czi";
maxPreviewImages = 3; // how many files to preview at once

biomassThresholdMethod = "Default";
biomassDarkBackground = false; // flip to true if the mask looks inverted
biomassGaussianSigma = 3;      // blurs out internal cell texture (organelles) before thresholding -- increase if mask still picks up internal detail
biomassMinParticleSize = 20;   // pixels^2 -- discards specks/dust smaller than this after masking; increase if dust still shows up, decrease if real small cells get removed

previewDir = getDirectory("Choose a folder with a few test .czi files (e.g. one condition folder, ideally from 2.1 or 3.1 so GFP signal is visible)");
allFiles = getFileList(previewDir);

// ---- If a channel window is a Z-stack, keep only the middle slice ----
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

count = 0;
for (f = 0; f < allFiles.length; f++) {
    if (count >= maxPreviewImages) break;
    if (!endsWith(allFiles[f], fileExtension)) continue;

    filePath = previewDir + allFiles[f];
    shortName = replace(allFiles[f], fileExtension, "");

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
        print("File " + allFiles[f] + " only has 1 channel -- skipping preview for this file.");
        continue;
    }

    if (c1Title == "" || c2Title == "") {
        print("Warning: could not find both C1 and C2 in " + allFiles[f] + " -- skipping.");
        continue;
    }

    // If either channel is a Z-stack, keep only the middle slice
    reduceToMiddleSlice(c1Title);
    reduceToMiddleSlice(c2Title);

    // Make raw display copies before thresholding destroys C1's pixel values visually
    selectWindow(c1Title);
    run("Duplicate...", "title=" + shortName + "_C1_raw duplicate");
    selectWindow(c2Title);
    run("Duplicate...", "title=" + shortName + "_C2_raw duplicate");

    // Segment C1 into the biomass mask:
    // 1. Blur out internal cell texture (organelles) so thresholding only
    //    catches the true cell/hyphae boundary, not internal detail.
    // 2. Threshold + convert to a binary mask.
    // 3. Fill Holes so cell interiors are solid, not hollow outlines.
    // 4. Analyze Particles with a minimum size filter to discard dust/specks.
    selectWindow(c1Title);
    run("Gaussian Blur...", "sigma=" + biomassGaussianSigma);
    if (biomassDarkBackground) {
        setAutoThreshold(biomassThresholdMethod + " dark");
    } else {
        setAutoThreshold(biomassThresholdMethod);
    }
    run("Convert to Mask");
    run("Fill Holes");
    run("Analyze Particles...", "size=" + biomassMinParticleSize + "-Infinity show=Masks");
    cleanMaskTitle = getTitle();
    run("Create Selection");

    if (selectionType() == -1) {
        print(allFiles[f] + ": no biomass detected with current settings.");
        biomassArea = 0;
    } else {
        biomassArea = getValue("Area");
        print(allFiles[f] + ": biomass area = " + biomassArea + " pixels");

        // Overlay the mask outline on both raw channels for visual comparison
        selectWindow(shortName + "_C1_raw");
        run("Restore Selection");
        run("Add Selection...");
        run("Select None");

        selectWindow(shortName + "_C2_raw");
        run("Restore Selection");
        run("Add Selection...");
        run("Select None");
    }

    close(); // closes cleanMaskTitle (the active window)
    selectWindow(c1Title);
    close(); // close the intermediate thresholded/masked version, keep the raw+overlay duplicates
    selectWindow(c2Title);
    close();

    // Arrange the two preview windows side by side for easy viewing
    selectWindow(shortName + "_C1_raw");
    setLocation(50, 100);
    selectWindow(shortName + "_C2_raw");
    setLocation(500, 100);

    count++;
}

print("Preview done. Check the '_C1_raw' and '_C2_raw' windows: the yellow outline");
print("shows the biomass mask boundary. It should trace around fungal material");
print("(cells, spores, or mycelium) -- not empty background.");
print("If it looks inverted, set biomassDarkBackground to the opposite value and re-run.");
