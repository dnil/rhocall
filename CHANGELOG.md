# ChangeLog

## unreleased
- HW and AZ are flags, set Number appropriately in VCF

## [0.5.1]
- fixed viz functionality for later versions of bcftools
- viz can now produce roh bed files

## [0.5]
- merged viz into the main rhocall framework
- viz can now produce wig files

## [0.4]
New flag --v14 to accept new bcftools roh files.

## [0.3]
--version

## [0.2]

Single chromosomes. Still some edge cases left for non-standard operation.
Classification pending. A test set has been added, but tests not yet written.
Looking at current test set performance, advise the use of bcftools v1.2 rather
than v.1.3 which seems to have broken something in their ROH functionality
(see bcftools #529).

## [0.1]

Initial release! Passes auto-install and performs intended functions on a
limited number of test samples.
Many edge cases can remain. Classification of AZ variant types is only
a placeholder right now. Future mode will include reading ped file,
noting inheritance and sex.
