import argparse
import glob
from pathlib import Path

import imageio

import cv2 as cv
import mitsuba

from tonemapper import *
from tqdm import tqdm

mitsuba.set_variant("scalar_rgb")


def read_steadyimg_mitsuba(filename: str) -> np.array:
    from mitsuba.core import Bitmap, Struct, float_dtype
    other = Bitmap(filename) \
        .convert(Bitmap.PixelFormat.RGBA, Struct.Type.Float32, srgb_gamma=False)
    img = np.array(other, copy=False)
    return img[:, :, :3]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Tool for Steady Images")
    # Options for reading steady and streak image
    parser.add_argument('-f', '--file', type=str, help="Steady image path",
                        default="image")
    args = parser.parse_args()

    img = read_steadyimg_mitsuba(args.file)
    img_tonemapped = tonemap(img,
                             normalize=True,
                             exposure=0,
                             offset=0,
                             tonemapper="GAMMA",
                             gamma=1.)
    imageio.imwrite(args.file.replace('.exr', '_tm.png'), (img_tonemapped * 255).astype(np.uint8))