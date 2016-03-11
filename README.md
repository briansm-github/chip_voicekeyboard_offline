# chip_voicekeyboard_online
Just some code for the 'offline' version of a USB voice keyboard
for NTC CHIP, as shown in video:
https://youtu.be/VgyPFYderjc

It is mostly bash shell script with some C to do the audio recording.

Note that the stock kernel does'nt have the necessary modules to support
USB HID, so I use the renzo kernel:
https://bbs.nextthing.co/t/compile-the-linux-kernel-for-chip-my-personal-howto/2669/40

    cd /tmp
    wget http://www.raspibo.org/renzo/chiplinux4.3.0rd235+.tgz
    cd /
    sudo tar zxf /tmp/chiplinux4.3.0rd235+.tgz

    ... and if you want to set this as your standard boot

    sudo cp /boot/vmlinuz-4.3.0rd235+ /boot/zImage

Also, need to disable the USB console as it conflicts otherwise

    cd /etc/modprobe.d
    echo "blacklist g_serial" > g_serial_blacklist.conf 

The entry point is init_offline, which can be run as an alternative init to systemd at boot-time, or run manually in a normally running CHIP.

The 'make_hid' script initializes the USB HID keyboard, creating
/dev/hidg0 device file to write to in order to generate keystrokes on the host.

The record/ code contains a very experimental recognizer using the FFT of LPC10 filter data as frame features. To train it, run something like
./go 26
and read the alphabet (a b c d etc) into the microphome. Training should only take a minute or so.
Then run the program with no arguments to test the recognition accuracy.

Some extra keystrokes are included in this version (6 more), to use them,
do:
./go 32
and after 'z' continue with 'space fullstop newline no on off'.

The code uses the Speex DSP library to remove as much noise as possible
from the recorded audio to improve quality.
