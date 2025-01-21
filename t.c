#include <stdio.h>
#include <stdlib.h>

int main() {
    // Command to change the wallpaper using gsettings for GNOME desktop
    const char* command = "gsettings set org.gnome.desktop.background picture-uri-dark 'file:///home/ian/Pictures/Screenshots/a.png'";

    // Execute the command
    int result = system(command);

    // Check if the command was executed successfully
    if (result == 0) {
        printf("Wallpaper changed successfully!\n");
    } else {
        printf("Failed to change the wallpaper.\n");
    }

    return 0;
}

