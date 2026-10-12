#[path = "charge_labels.rs"]
mod charge_labels;

/// Match the charge left after CCTBX's scattering-label normalization.
/// Unsupported ionic scattering labels become their neutral element label.
pub(super) fn normalized(element: &str, charge: i32) -> i32 {
    if charge == 0 {
        return 0;
    }
    let label = format!(
        "{}{}{}",
        element.trim().to_ascii_uppercase(),
        charge.abs(),
        if charge > 0 { '+' } else { '-' }
    );
    if charge_labels::IONIC_LABELS
        .binary_search(&label.as_str())
        .is_ok()
    {
        charge
    } else {
        0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn only_supported_ionic_forms_keep_their_charge() {
        assert_eq!(normalized("H", 1), 0);
        assert_eq!(normalized("D", 1), 0);
        assert_eq!(normalized("H", -1), -1);
        assert_eq!(normalized("O", -1), -1);
        assert_eq!(normalized("O", -2), -2);
        assert_eq!(normalized("N", 1), 0);
        assert_eq!(normalized("Fe", 2), 2);
        assert_eq!(normalized("Fe", 9), 0);
    }
}
