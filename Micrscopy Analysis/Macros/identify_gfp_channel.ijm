// Diagnostic: open a single .czi file, split its channels, and display
// each one individually with its index labeled in the title bar, so we
// can identify which channel number corresponds to GFP.

path = File.openDialog("Choose one representative .czi file (ideally from a 2.1 or 3.1 folder, since it should show visible GFP signal)");

run("Bio-Formats Importer", "open=[" + path + "] color_mode=Default view=Hyperstack stack_order=XYCZT");

title = getTitle();
print("Opened: " + title);
getDimensions(width, height, channels, slices, frames);
print("Width=" + width + " Height=" + height + " Channels=" + channels + " Slices=" + slices + " Frames=" + frames);

if (channels > 1) {
    run("Split Channels");
    print("Channels split into separate windows -- inspect each one below.");
    print("Look at the title of each window (e.g. 'C1-...', 'C2-...') and note which one shows clear GFP/fluorescent signal.");
} else {
    print("This file only has 1 channel -- no splitting needed.");
}
