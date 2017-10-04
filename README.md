## MVCT Retrospective Analysis Tool

by Mark Geurts <mark.w.geurts@gmail.com>
<br>Copyright &copy; 2017, University of Wisconsin Board of Regents

The MVCT Retrospective Analysis Tool is a command-line application written in MATLAB&reg; that searches a folder for [TomoTherapy&reg;](http://www.accuray.com) patient archives, finds performed MVCT scans, and writes information about each scan to a .csv file for evaluation. In addition to extracting information about each scan, the tool can optionally also re-register the images using MATLAB or [plastimatch<sup>&copy;</sup>](http://plastimatch.org) 3D/6D rigid registration and report differences in image similarity.

TomoTherapy is a registered trademark of Accuray Incorporated. MATLAB is a registered trademark of MathWorks Inc. Plastimatch is Copyrighted by The General Hospital Corporation and is available at http://plastimatch.org.

## Contents

* [Installation and Use](README.md#installation-and-use)
* [Configuration File](README.md#configuration-file)
* [Report File](README.md#report-file)
* [Compatibility and Requirements](README.md#compatibility-and-requirements)
* [Troubleshooting](README.md#troubleshooting)
* [License](README.md#license) 

## Installation and Use

To install the MVCT Retrospective Analysis Tool, navigate the target folder and run `git clone --recursive https://github.com/mwgeurts/mvct_retro` from a terminal compatible with git. This will download all files from this repository as well as included submodule repositories. Once downloaded, the application can be configured by modifying the `config.txt` configuration file and run from MATLAB by executing `ScanArchives`.

Because this application is command-line only, a graphical environment is not required. In fact, it is recommended that MATLAB be started with display features disabled, such as `matlab -noFigureWindows` on Windows.

Finally, when executing this tool on large directories of archives, the tool can be stopped and re-started without re-analyzing or duplicating existing data. The tool will read in the current results .csv file and skip any scan where the same scan UID and tool version exists.

## Configuration File

The configuration file `config.txt` contains various settings that can be set prior to tool execution. The following table describes each input and its default setting (if applicable).

| Configuration Setting | Default Value | Description |
|-----------------------|---------------|-------------|
| ARCHIVE_PATH |  | This is the directory that `ScanArchives()` will search for patient archives. All subfolders are scanned recursively. The path can defined relative to the current MATLAB path or as an absolute location. |
| RESULTS_CSV | `../results.csv` | This is the path and filename of the results .csv file. If this file does not exist, the tool will automatically create it. |
| ANON_RESULTS | `0` | Flag indicating whether to anonymize the results.csv file `1` or not `0`. When anonymized, the archive folder path will not be written to the results .csv file (most archive folders contain the patient name). |
| REGISTRATION_METHOD | `MATLAB` | If this configuration option is provided, the tool will re-register each MVCT to the planning CT and record the registration adjustments. The value can be `MATLAB` or `PLASTIMATCH`. See [fusion_tools](https://github.com/mwgeurts/fusion_tools) for more information. |
| REGISTRATION_METRIC | `MI` | If `REGISTRATION_METHOD` is provided, the metric used to re-register the images. Can be `MSE`, `MI`, or (plastimatch only) `GM`. |
| REGISTRATION_LEVELS | `3` | If `REGISTRATION_METHOD` is provided, the number of downsampling levels to run during registration. |
| REGISTRATION_ITER | `50` | If `REGISTRATION_METHOD` is provided, the maximum number of iterations to run at each level. |
| REGISTER_BONE | '0' | If `REGISTRATION_METHOD` is provided, a flag indicating whether to re-register each MVCT scan using a bone-only algorithm `1` or to use the full image `0`. |
| ALLOW_ROTATIONS | `1` | If `REGISTRATION_METHOD` is provided, a flag indicating whether to include rotations `1` or only translations `0`. |
| SIMILARITY_METRIC | `SSI` | If this configuration option is provided, the tool will calculate the image similarity of the user-registration merged MVCT. If provided with a `REGISTRATION_METHOD` option, the re-registered image similarity will also be reported. Values can be `MSE` or `SSI` |

## Report File

The RESULTS_CSV file contains the following columns. The registration and
similarity columns will only be filled out if set in `config.txt`.

| Column | Value |
|--------|-------|
| 1 | Full path to patient archive _patient.xml.  However, if the config option ANON_RESULTS is set to 1, will be empty. |
| 2 |  SHA1 signature of _patient.xml file |
| 3 |  MVCT Timestamp (MATLAB datenum) |
| 4 |  MVCT UID |
| 5 |  Plan Name |
| 6 |  Scan Length (cm) |
| 7 |  Time from scan to treatment (minutes) |
| 8 |  Number of MVCT scans performed on same day |
| 9 |  User Registered Pitch (degrees) |
| 10 |  User Registered Yaw (degrees) |
| 11 |  User Registered Roll (degrees) |
| 12 | User Registered X Translation (cm) |
| 13 | User Registered Y Translation (cm) |
| 14 | User Registered Z Translation (cm) |
| 15 | Tool Version |
| 16 | Registration Method, if provided |
| 17 | Re-Registered Pitch (degrees) |
| 18 | Re-Registered Yaw (degrees) |
| 19 | Re-Registered Roll (degrees) |
| 20 | Re-Registered X Translation (cm) |
| 21 | Re-Registered Y Translation (cm) |
| 22 | Re-Registered Z Translation (cm) |
| 23 | Similarity Metric, if provided |
| 24 | User Registration Similarity |
| 25 | Re-Registration Similarity |

## Compatibility and Requirements

The MVCT Retrospective Analysis Tool has been validated for MATLAB versions 9.1 through 9.2 on Macintosh OSX 10.10 (Yosemite), Windows 7, and Ubuntu 14.04. These tools optionally use the Image Processing Toolbox MATLAB function `imregtform()` when performing MATLAB-based rigid registration and the Parallel Processing Toolbox `interp3()` when performing GPU-based interpolation. If this toolbox or a compatible GPU device is not found, the `MergeImages()` function will automatically revert to CPU interpolation. Plastimatch registration does not require any additional toolboxes.

## Troubleshooting

This application records key input parameters and results to a log.txt file using the `Event()` function. The log is the most important route to troubleshooting errors encountered by this software. **Warning! the log file may contain protected health information.** The author can also be contacted using the information above. Refer to the license file for a full description of the limitations on liability when using or this software or its components.

## License

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.

