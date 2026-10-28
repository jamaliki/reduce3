//! Where the tests find chem_data (`REDUCE3_CHEM_DATA`, see
//! [`reduce3::MonLib::locate`]). A test that needs it is skipped when it is
//! absent, unless `REDUCE3_REQUIRE_CHEM_DATA` is set: then the test fails, so
//! a run that is meant to cover them (CI) cannot pass by skipping them.

use std::path::PathBuf;

/// Set to anything but empty or `0` to fail, not skip, without chem_data.
pub const REQUIRE: &str = "REDUCE3_REQUIRE_CHEM_DATA";

pub fn required() -> bool {
    std::env::var_os(REQUIRE).is_some_and(|v| !v.is_empty() && v != "0")
}

/// The chem_data root, or None to skip the test.
pub fn locate() -> Option<PathBuf> {
    let root = reduce3::MonLib::locate(None);
    if root.is_none() {
        assert!(!required(), "chem_data not found and {REQUIRE} is set: set REDUCE3_CHEM_DATA");
        eprintln!("skipped: chem_data not found (set REDUCE3_CHEM_DATA)");
    }
    root
}
