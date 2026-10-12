# Archiving corrected versions

The concept DOI [10.5281/zenodo.20746987](https://doi.org/10.5281/zenodo.20746987)
identifies the evolving project and resolves to its latest published Zenodo version.
The version DOI [10.5281/zenodo.23055573](https://doi.org/10.5281/zenodo.23055573)
identifies the published v0.4.0 archive. Use a version DOI when citing the exact
software or manuscript used in an experiment.

A local edit or GitHub commit does not update files already archived on Zenodo.
Substantive corrections to proofs, reproduction scripts and manuscripts should be
archived as a new linked version, with a new version DOI. Keep the concept DOI;
do not guess the identifier Zenodo will assign to the new version. The existing
version/date fields in CITATION.cff describe the last published release until a
new release is prepared.

Zenodo permits metadata edits without changing the DOI. Its current file policy
also permits minor file corrections within 30 days of publication, with that edit
draft published within 45 days of the original publication. This exception is
separate from versioning and does not make GitHub commits update an archive.

Sources: [DOI versioning](https://zenodo.org/help/versioning),
[manage versions](https://help.zenodo.org/docs/deposit/manage-versions/),
[manage files](https://help.zenodo.org/docs/deposit/manage-files/).

Before publication, review the exact repository patch, manuscript PDFs, release
version, release notes and archive destination. Publish only after the owner
explicitly approves that concrete content and destination. A GitHub release may
trigger archive automation if the repository is connected to Zenodo; verify that
connection and the resulting archived files rather than assuming it exists.
