//! Color definition and conversion — analogous to OCCT `Quantity_Color` / `Quantity_ColorRGBA`.
//!
//! Supports color spaces: linear RGB, sRGB, HLS, CIELab, CIELch.
//! Internal storage is linear RGB `[f32; 3]` (matching OCCT's `NCollection_Vec3<float>`).

use std::fmt;

// ============================================================================
// TypeOfColor (OCCT Quantity_TypeOfColor)
// ============================================================================

/// Color definition systems.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TypeOfColor {
    /// Normalized linear RGB, each component in [0, 1].
    RGB,
    /// Non-linear gamma-shifted sRGB, each component in [0, 1].
    sRGB,
    /// Hue (degrees [0, 360], -1 for gray) / Lightness [0, 1] / Saturation [0, 1].
    HLS,
    /// CIE L*a*b* with D65 white point: L [0, 100], a/b approx [-110, 100].
    CIELab,
    /// CIE L*c*h*: L [0, 100], c [0, ~135], h [0, 360].
    CIELch,
}

// ============================================================================
// NameOfColor (OCCT Quantity_NameOfColor)
// ============================================================================

/// Named colors (X11 specification based, matching OCCT).
/// Only the most common subset is exposed as enum variants.
#[allow(missing_docs)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum NameOfColor {
    Black,
    Blue,
    Brown,
    Cyan,
    DarkGreen,
    DarkOrange,
    DarkRed,
    Gold,
    Gray,
    Green,
    Ivory,
    Khaki,
    Lavender,
    LimeGreen,
    Magenta,
    Maroon,
    NavyBlue,
    Olive,
    Orange,
    Orchid,
    Pink,
    Plum,
    Purple,
    Red,
    Salmon,
    SeaGreen,
    Sienna,
    Silver,
    SkyBlue,
    SlateGray,
    Snow,
    SteelBlue,
    Tan,
    Teal,
    Tomato,
    Turquoise,
    Violet,
    Wheat,
    White,
    Yellow,
    AliceBlue,
    AntiqueWhite,
    Aquamarine,
    Azure,
    Beige,
    Bisque,
    BlanchedAlmond,
    BlueViolet,
    BurlyWood,
    CadetBlue,
    Chartreuse,
    Chocolate,
    Coral,
    CornflowerBlue,
    Cornsilk,
    DarkGoldenrod,
    DarkKhaki,
    DarkOliveGreen,
    DarkOrchid,
    DarkSalmon,
    DarkSeaGreen,
    DarkSlateBlue,
    DarkSlateGray,
    DarkTurquoise,
    DarkViolet,
    DeepPink,
    DeepSkyBlue,
    DimGray,
    DodgerBlue,
    Firebrick,
    FloralWhite,
    ForestGreen,
    Fuchsia,
    Gainsboro,
    GhostWhite,
    Goldenrod,
    GreenYellow,
    Honeydew,
    HotPink,
    IndianRed,
    Indigo,
    LavenderBlush,
    LawnGreen,
    LemonChiffon,
    LightBlue,
    LightCoral,
    LightCyan,
    LightGoldenrod,
    LightGoldenrodYellow,
    LightGray,
    LightPink,
    LightSalmon,
    LightSeaGreen,
    LightSkyBlue,
    LightSlateBlue,
    LightSlateGray,
    LightSteelBlue,
    LightYellow,
    Linen,
    MediumAquamarine,
    MediumBlue,
    MediumOrchid,
    MediumPurple,
    MediumSeaGreen,
    MediumSlateBlue,
    MediumSpringGreen,
    MediumTurquoise,
    MediumVioletRed,
    MidnightBlue,
    MintCream,
    MistyRose,
    Moccasin,
    NavajoWhite,
    OldLace,
    OliveDrab,
    OrangeRed,
    PaleGoldenrod,
    PaleGreen,
    PaleTurquoise,
    PaleVioletRed,
    PapayaWhip,
    PeachPuff,
    Peru,
    PowderBlue,
    RosyBrown,
    RoyalBlue,
    SaddleBrown,
    SandyBrown,
    Seashell,
    SlateBlue,
    SpringGreen,
    Thistle,
    VioletRed,
    WhiteSmoke,
    YellowGreen,
}

impl NameOfColor {
    /// OCCT: `Quantity_Color::StringName(theColor)`.
    pub fn string_name(self) -> &'static str {
        match self {
            NameOfColor::Black => "BLACK",
            NameOfColor::Blue => "BLUE",
            NameOfColor::Brown => "BROWN",
            NameOfColor::Cyan => "CYAN",
            NameOfColor::DarkGreen => "DARK_GREEN",
            NameOfColor::DarkOrange => "DARK_ORANGE",
            NameOfColor::DarkRed => "DARK_RED",
            NameOfColor::Gold => "GOLD",
            NameOfColor::Gray => "GRAY",
            NameOfColor::Green => "GREEN",
            NameOfColor::Ivory => "IVORY",
            NameOfColor::Khaki => "KHAKI",
            NameOfColor::Lavender => "LAVENDER",
            NameOfColor::LimeGreen => "LIME_GREEN",
            NameOfColor::Magenta => "MAGENTA",
            NameOfColor::Maroon => "MAROON",
            NameOfColor::NavyBlue => "NAVY_BLUE",
            NameOfColor::Olive => "OLIVE",
            NameOfColor::Orange => "ORANGE",
            NameOfColor::Orchid => "ORCHID",
            NameOfColor::Pink => "PINK",
            NameOfColor::Plum => "PLUM",
            NameOfColor::Purple => "PURPLE",
            NameOfColor::Red => "RED",
            NameOfColor::Salmon => "SALMON",
            NameOfColor::SeaGreen => "SEA_GREEN",
            NameOfColor::Sienna => "SIENNA",
            NameOfColor::Silver => "SILVER",
            NameOfColor::SkyBlue => "SKY_BLUE",
            NameOfColor::SlateGray => "SLATE_GRAY",
            NameOfColor::Snow => "SNOW",
            NameOfColor::SteelBlue => "STEEL_BLUE",
            NameOfColor::Tan => "TAN",
            NameOfColor::Teal => "TEAL",
            NameOfColor::Tomato => "TOMATO",
            NameOfColor::Turquoise => "TURQUOISE",
            NameOfColor::Violet => "VIOLET",
            NameOfColor::Wheat => "WHEAT",
            NameOfColor::White => "WHITE",
            NameOfColor::Yellow => "YELLOW",
            NameOfColor::AliceBlue => "ALICE_BLUE",
            NameOfColor::AntiqueWhite => "ANTIQUE_WHITE",
            NameOfColor::Aquamarine => "AQUAMARINE",
            NameOfColor::Azure => "AZURE",
            NameOfColor::Beige => "BEIGE",
            NameOfColor::Bisque => "BISQUE",
            NameOfColor::BlanchedAlmond => "BLANCHED_ALMOND",
            NameOfColor::BlueViolet => "BLUE_VIOLET",
            NameOfColor::BurlyWood => "BURLY_WOOD",
            NameOfColor::CadetBlue => "CADET_BLUE",
            NameOfColor::Chartreuse => "CHARTREUSE",
            NameOfColor::Chocolate => "CHOCOLATE",
            NameOfColor::Coral => "CORAL",
            NameOfColor::CornflowerBlue => "CORNFLOWER_BLUE",
            NameOfColor::Cornsilk => "CORNSILK",
            NameOfColor::DarkGoldenrod => "DARK_GOLDENROD",
            NameOfColor::DarkKhaki => "DARK_KHAKI",
            NameOfColor::DarkOliveGreen => "DARK_OLIVE_GREEN",
            NameOfColor::DarkOrchid => "DARK_ORCHID",
            NameOfColor::DarkSalmon => "DARK_SALMON",
            NameOfColor::DarkSeaGreen => "DARK_SEA_GREEN",
            NameOfColor::DarkSlateBlue => "DARK_SLATE_BLUE",
            NameOfColor::DarkSlateGray => "DARK_SLATE_GRAY",
            NameOfColor::DarkTurquoise => "DARK_TURQUOISE",
            NameOfColor::DarkViolet => "DARK_VIOLET",
            NameOfColor::DeepPink => "DEEP_PINK",
            NameOfColor::DeepSkyBlue => "DEEP_SKY_BLUE",
            NameOfColor::DimGray => "DIM_GRAY",
            NameOfColor::DodgerBlue => "DODGER_BLUE",
            NameOfColor::Firebrick => "FIREBRICK",
            NameOfColor::FloralWhite => "FLORAL_WHITE",
            NameOfColor::ForestGreen => "FOREST_GREEN",
            NameOfColor::Fuchsia => "FUCHSIA",
            NameOfColor::Gainsboro => "GAINSBORO",
            NameOfColor::GhostWhite => "GHOST_WHITE",
            NameOfColor::Goldenrod => "GOLDENROD",
            NameOfColor::GreenYellow => "GREEN_YELLOW",
            NameOfColor::Honeydew => "HONEYDEW",
            NameOfColor::HotPink => "HOT_PINK",
            NameOfColor::IndianRed => "INDIAN_RED",
            NameOfColor::Indigo => "INDIGO",
            NameOfColor::LavenderBlush => "LAVENDER_BLUSH",
            NameOfColor::LawnGreen => "LAWN_GREEN",
            NameOfColor::LemonChiffon => "LEMON_CHIFFON",
            NameOfColor::LightBlue => "LIGHT_BLUE",
            NameOfColor::LightCoral => "LIGHT_CORAL",
            NameOfColor::LightCyan => "LIGHT_CYAN",
            NameOfColor::LightGoldenrod => "LIGHT_GOLDENROD",
            NameOfColor::LightGoldenrodYellow => "LIGHT_GOLDENROD_YELLOW",
            NameOfColor::LightGray => "LIGHT_GRAY",
            NameOfColor::LightPink => "LIGHT_PINK",
            NameOfColor::LightSalmon => "LIGHT_SALMON",
            NameOfColor::LightSeaGreen => "LIGHT_SEA_GREEN",
            NameOfColor::LightSkyBlue => "LIGHT_SKY_BLUE",
            NameOfColor::LightSlateBlue => "LIGHT_SLATE_BLUE",
            NameOfColor::LightSlateGray => "LIGHT_SLATE_GRAY",
            NameOfColor::LightSteelBlue => "LIGHT_STEEL_BLUE",
            NameOfColor::LightYellow => "LIGHT_YELLOW",
            NameOfColor::Linen => "LINEN",
            NameOfColor::MediumAquamarine => "MEDIUM_AQUAMARINE",
            NameOfColor::MediumBlue => "MEDIUM_BLUE",
            NameOfColor::MediumOrchid => "MEDIUM_ORCHID",
            NameOfColor::MediumPurple => "MEDIUM_PURPLE",
            NameOfColor::MediumSeaGreen => "MEDIUM_SEA_GREEN",
            NameOfColor::MediumSlateBlue => "MEDIUM_SLATE_BLUE",
            NameOfColor::MediumSpringGreen => "MEDIUM_SPRING_GREEN",
            NameOfColor::MediumTurquoise => "MEDIUM_TURQUOISE",
            NameOfColor::MediumVioletRed => "MEDIUM_VIOLET_RED",
            NameOfColor::MidnightBlue => "MIDNIGHT_BLUE",
            NameOfColor::MintCream => "MINT_CREAM",
            NameOfColor::MistyRose => "MISTY_ROSE",
            NameOfColor::Moccasin => "MOCCASIN",
            NameOfColor::NavajoWhite => "NAVAJO_WHITE",
            NameOfColor::OldLace => "OLD_LACE",
            NameOfColor::OliveDrab => "OLIVE_DRAB",
            NameOfColor::OrangeRed => "ORANGE_RED",
            NameOfColor::PaleGoldenrod => "PALE_GOLDENROD",
            NameOfColor::PaleGreen => "PALE_GREEN",
            NameOfColor::PaleTurquoise => "PALE_TURQUOISE",
            NameOfColor::PaleVioletRed => "PALE_VIOLET_RED",
            NameOfColor::PapayaWhip => "PAPAYA_WHIP",
            NameOfColor::PeachPuff => "PEACH_PUFF",
            NameOfColor::Peru => "PERU",
            NameOfColor::PowderBlue => "POWDER_BLUE",
            NameOfColor::RosyBrown => "ROSY_BROWN",
            NameOfColor::RoyalBlue => "ROYAL_BLUE",
            NameOfColor::SaddleBrown => "SADDLE_BROWN",
            NameOfColor::SandyBrown => "SANDY_BROWN",
            NameOfColor::Seashell => "SEASHELL",
            NameOfColor::SlateBlue => "SLATE_BLUE",
            NameOfColor::SpringGreen => "SPRING_GREEN",
            NameOfColor::Thistle => "THISTLE",
            NameOfColor::VioletRed => "VIOLET_RED",
            NameOfColor::WhiteSmoke => "WHITE_SMOKE",
            NameOfColor::YellowGreen => "YELLOW_GREEN",
        }
    }
}

// ============================================================================
// Named color RGB values
// ============================================================================

/// Lookup the linear RGB values for a named color.
fn named_color_rgb(name: NameOfColor) -> [f32; 3] {
    match name {
        NameOfColor::Black => [0.0, 0.0, 0.0],
        NameOfColor::Blue => [0.0, 0.0, 1.0],
        NameOfColor::Brown => [0.6471, 0.1647, 0.1647],
        NameOfColor::Cyan => [0.0, 1.0, 1.0],
        NameOfColor::DarkGreen => [0.0, 0.3922, 0.0],
        NameOfColor::DarkOrange => [1.0, 0.5490, 0.0],
        NameOfColor::DarkRed => [0.5451, 0.0, 0.0],
        NameOfColor::Gold => [1.0, 0.8431, 0.0],
        NameOfColor::Gray => [0.5020, 0.5020, 0.5020],
        NameOfColor::Green => [0.0, 0.5020, 0.0],
        NameOfColor::Ivory => [1.0, 1.0, 0.9412],
        NameOfColor::Khaki => [0.9412, 0.9020, 0.5490],
        NameOfColor::Lavender => [0.9020, 0.9020, 0.9804],
        NameOfColor::LimeGreen => [0.1961, 0.8039, 0.1961],
        NameOfColor::Magenta => [1.0, 0.0, 1.0],
        NameOfColor::Maroon => [0.5020, 0.0, 0.0],
        NameOfColor::NavyBlue => [0.0, 0.0, 0.5020],
        NameOfColor::Olive => [0.5020, 0.5020, 0.0],
        NameOfColor::Orange => [1.0, 0.6471, 0.0],
        NameOfColor::Orchid => [0.8549, 0.4392, 0.8392],
        NameOfColor::Pink => [1.0, 0.7529, 0.7961],
        NameOfColor::Plum => [0.8667, 0.6275, 0.8667],
        NameOfColor::Purple => [0.5020, 0.0, 0.5020],
        NameOfColor::Red => [1.0, 0.0, 0.0],
        NameOfColor::Salmon => [0.9804, 0.5020, 0.4471],
        NameOfColor::SeaGreen => [0.1804, 0.5451, 0.3412],
        NameOfColor::Sienna => [0.6275, 0.3216, 0.1765],
        NameOfColor::Silver => [0.7529, 0.7529, 0.7529],
        NameOfColor::SkyBlue => [0.5294, 0.8078, 0.9216],
        NameOfColor::SlateGray => [0.4392, 0.5020, 0.5647],
        NameOfColor::Snow => [1.0, 0.9804, 0.9804],
        NameOfColor::SteelBlue => [0.2745, 0.5098, 0.7059],
        NameOfColor::Tan => [0.8235, 0.7059, 0.5490],
        NameOfColor::Teal => [0.0, 0.5020, 0.5020],
        NameOfColor::Tomato => [1.0, 0.3882, 0.2784],
        NameOfColor::Turquoise => [0.2510, 0.8784, 0.8157],
        NameOfColor::Violet => [0.9333, 0.5098, 0.9333],
        NameOfColor::Wheat => [0.9608, 0.8706, 0.7020],
        NameOfColor::White => [1.0, 1.0, 1.0],
        NameOfColor::Yellow => [1.0, 1.0, 0.0],
        NameOfColor::AliceBlue => [0.9412, 0.9725, 1.0],
        NameOfColor::AntiqueWhite => [0.9804, 0.9216, 0.8431],
        NameOfColor::Aquamarine => [0.4980, 1.0, 0.8314],
        NameOfColor::Azure => [0.9412, 1.0, 1.0],
        NameOfColor::Beige => [0.9608, 0.9608, 0.8627],
        NameOfColor::Bisque => [1.0, 0.8941, 0.7686],
        NameOfColor::BlanchedAlmond => [1.0, 0.9216, 0.8039],
        NameOfColor::BlueViolet => [0.5412, 0.1686, 0.8863],
        NameOfColor::BurlyWood => [0.8706, 0.7216, 0.5294],
        NameOfColor::CadetBlue => [0.3725, 0.6196, 0.6275],
        NameOfColor::Chartreuse => [0.4980, 1.0, 0.0],
        NameOfColor::Chocolate => [0.8235, 0.4118, 0.1176],
        NameOfColor::Coral => [1.0, 0.4980, 0.3137],
        NameOfColor::CornflowerBlue => [0.3922, 0.5843, 0.9294],
        NameOfColor::Cornsilk => [1.0, 0.9725, 0.8627],
        NameOfColor::DarkGoldenrod => [0.7216, 0.5255, 0.0431],
        NameOfColor::DarkKhaki => [0.7412, 0.7176, 0.4196],
        NameOfColor::DarkOliveGreen => [0.3333, 0.4196, 0.1843],
        NameOfColor::DarkOrchid => [0.6000, 0.1961, 0.8000],
        NameOfColor::DarkSalmon => [0.9137, 0.5882, 0.4784],
        NameOfColor::DarkSeaGreen => [0.5608, 0.7373, 0.5608],
        NameOfColor::DarkSlateBlue => [0.2824, 0.2392, 0.5451],
        NameOfColor::DarkSlateGray => [0.1843, 0.3098, 0.3098],
        NameOfColor::DarkTurquoise => [0.0, 0.8078, 0.8196],
        NameOfColor::DarkViolet => [0.5804, 0.0, 0.8275],
        NameOfColor::DeepPink => [1.0, 0.0784, 0.5765],
        NameOfColor::DeepSkyBlue => [0.0, 0.7490, 1.0],
        NameOfColor::DimGray => [0.4118, 0.4118, 0.4118],
        NameOfColor::DodgerBlue => [0.1176, 0.5647, 1.0],
        NameOfColor::Firebrick => [0.6980, 0.1333, 0.1333],
        NameOfColor::FloralWhite => [1.0, 0.9804, 0.9412],
        NameOfColor::ForestGreen => [0.1333, 0.5451, 0.1333],
        NameOfColor::Fuchsia => [1.0, 0.0, 1.0],
        NameOfColor::Gainsboro => [0.8627, 0.8627, 0.8627],
        NameOfColor::GhostWhite => [0.9725, 0.9725, 1.0],
        NameOfColor::Goldenrod => [0.8549, 0.6471, 0.1255],
        NameOfColor::GreenYellow => [0.6784, 1.0, 0.1843],
        NameOfColor::Honeydew => [0.9412, 1.0, 0.9412],
        NameOfColor::HotPink => [1.0, 0.4118, 0.7059],
        NameOfColor::IndianRed => [0.8039, 0.3608, 0.3608],
        NameOfColor::Indigo => [0.2941, 0.0, 0.5098],
        NameOfColor::LavenderBlush => [1.0, 0.9412, 0.9608],
        NameOfColor::LawnGreen => [0.4863, 0.9882, 0.0],
        NameOfColor::LemonChiffon => [1.0, 0.9804, 0.8039],
        NameOfColor::LightBlue => [0.6784, 0.8471, 0.9020],
        NameOfColor::LightCoral => [0.9412, 0.5020, 0.5020],
        NameOfColor::LightCyan => [0.8784, 1.0, 1.0],
        NameOfColor::LightGoldenrod => [0.9333, 0.8667, 0.5098],
        NameOfColor::LightGoldenrodYellow => [0.9804, 0.9804, 0.8235],
        NameOfColor::LightGray => [0.8275, 0.8275, 0.8275],
        NameOfColor::LightPink => [1.0, 0.7137, 0.7569],
        NameOfColor::LightSalmon => [1.0, 0.6275, 0.4784],
        NameOfColor::LightSeaGreen => [0.1255, 0.6980, 0.6667],
        NameOfColor::LightSkyBlue => [0.5294, 0.8078, 0.9804],
        NameOfColor::LightSlateBlue => [0.5176, 0.4392, 1.0],
        NameOfColor::LightSlateGray => [0.4667, 0.5333, 0.6000],
        NameOfColor::LightSteelBlue => [0.6902, 0.7686, 0.8706],
        NameOfColor::LightYellow => [1.0, 1.0, 0.8784],
        NameOfColor::Linen => [0.9804, 0.9412, 0.9020],
        NameOfColor::MediumAquamarine => [0.4000, 0.8039, 0.6667],
        NameOfColor::MediumBlue => [0.0, 0.0, 0.8039],
        NameOfColor::MediumOrchid => [0.7294, 0.3333, 0.8275],
        NameOfColor::MediumPurple => [0.5765, 0.4392, 0.8588],
        NameOfColor::MediumSeaGreen => [0.2353, 0.7020, 0.4431],
        NameOfColor::MediumSlateBlue => [0.4824, 0.4078, 0.9333],
        NameOfColor::MediumSpringGreen => [0.0, 0.9804, 0.6039],
        NameOfColor::MediumTurquoise => [0.2824, 0.8196, 0.8000],
        NameOfColor::MediumVioletRed => [0.7804, 0.0824, 0.5216],
        NameOfColor::MidnightBlue => [0.0980, 0.0980, 0.4392],
        NameOfColor::MintCream => [0.9608, 1.0, 0.9804],
        NameOfColor::MistyRose => [1.0, 0.8941, 0.8824],
        NameOfColor::Moccasin => [1.0, 0.8941, 0.7098],
        NameOfColor::NavajoWhite => [1.0, 0.8706, 0.6784],
        NameOfColor::OldLace => [0.9922, 0.9608, 0.9020],
        NameOfColor::OliveDrab => [0.4196, 0.5569, 0.1373],
        NameOfColor::OrangeRed => [1.0, 0.2706, 0.0],
        NameOfColor::PaleGoldenrod => [0.9333, 0.9098, 0.6667],
        NameOfColor::PaleGreen => [0.5961, 0.9843, 0.5961],
        NameOfColor::PaleTurquoise => [0.6863, 0.9333, 0.9333],
        NameOfColor::PaleVioletRed => [0.8588, 0.4392, 0.5765],
        NameOfColor::PapayaWhip => [1.0, 0.9373, 0.8353],
        NameOfColor::PeachPuff => [1.0, 0.8549, 0.7255],
        NameOfColor::Peru => [0.8039, 0.5216, 0.2471],
        NameOfColor::PowderBlue => [0.6902, 0.8784, 0.9020],
        NameOfColor::RosyBrown => [0.7373, 0.5608, 0.5608],
        NameOfColor::RoyalBlue => [0.2549, 0.4118, 0.8824],
        NameOfColor::SaddleBrown => [0.5451, 0.2706, 0.0745],
        NameOfColor::SandyBrown => [0.9569, 0.6431, 0.3765],
        NameOfColor::Seashell => [1.0, 0.9608, 0.9333],
        NameOfColor::SlateBlue => [0.4157, 0.3529, 0.8039],
        NameOfColor::SpringGreen => [0.0, 1.0, 0.4980],
        NameOfColor::Thistle => [0.8471, 0.7490, 0.8471],
        NameOfColor::VioletRed => [0.8157, 0.1255, 0.5647],
        NameOfColor::WhiteSmoke => [0.9608, 0.9608, 0.9608],
        NameOfColor::YellowGreen => [0.6039, 0.8039, 0.1961],
    }
}

// ============================================================================
// Color space conversion helpers
// ============================================================================

/// Linear RGB → sRGB (gamma correction per OpenGL specs).
fn linear_to_srgb(c: f32) -> f32 {
    if c <= 0.0031308 {
        c * 12.92
    } else {
        c.powf(1.0 / 2.4) * 1.055 - 0.055
    }
}

/// sRGB → linear RGB.
fn srgb_to_linear(c: f32) -> f32 {
    if c <= 0.04045 {
        c / 12.92
    } else {
        ((c + 0.055) / 1.055).powf(2.4)
    }
}

fn linear_to_srgb_v3(c: [f32; 3]) -> [f32; 3] {
    [linear_to_srgb(c[0]), linear_to_srgb(c[1]), linear_to_srgb(c[2])]
}

fn srgb_to_linear_v3(c: [f32; 3]) -> [f32; 3] {
    [srgb_to_linear(c[0]), srgb_to_linear(c[1]), srgb_to_linear(c[2])]
}

/// sRGB → HLS (OCCT Quantity_Color::Convert_sRGB_To_HLS).
fn srgb_to_hls(rgb: [f32; 3]) -> [f32; 3] {
    let r = rgb[0] as f64;
    let g = rgb[1] as f64;
    let b = rgb[2] as f64;

    let mn = r.min(g).min(b);
    let mx = r.max(g).max(b);
    let l = (mn + mx) / 2.0;

    let delta = mx - mn;
    if delta < 1e-12 {
        return [-1.0, l as f32, 0.0]; // achromatic
    }

    let s = if l <= 0.5 {
        delta / (mx + mn)
    } else {
        delta / (2.0 - mx - mn)
    };

    let mut h = if (mx - r).abs() < 1e-12 {
        (g - b) / delta
    } else if (mx - g).abs() < 1e-12 {
        2.0 + (b - r) / delta
    } else {
        4.0 + (r - g) / delta
    };

    h *= 60.0;
    if h < 0.0 {
        h += 360.0;
    }

    [h as f32, l as f32, s as f32]
}

/// HLS → sRGB (OCCT Quantity_Color::Convert_HLS_To_sRGB).
fn hls_to_srgb(hls: [f32; 3]) -> [f32; 3] {
    let h = hls[0] as f64;
    let l = hls[1] as f64;
    let s = hls[2] as f64;

    if s < 1e-12 {
        let v = l.max(0.0).min(1.0);
        return [v as f32, v as f32, v as f32];
    }

    let m2 = if l <= 0.5 {
        l * (1.0 + s)
    } else {
        l + s - l * s
    };
    let m1 = 2.0 * l - m2;

    let hue = h / 360.0;
    let r = hue_value(m1, m2, hue + 1.0 / 3.0);
    let g = hue_value(m1, m2, hue);
    let b = hue_value(m1, m2, hue - 1.0 / 3.0);

    [r as f32, g as f32, b as f32]
}

fn hue_value(m1: f64, m2: f64, mut h: f64) -> f64 {
    if h < 0.0 {
        h += 1.0;
    }
    if h > 1.0 {
        h -= 1.0;
    }
    if h < 1.0 / 6.0 {
        m1 + (m2 - m1) * h * 6.0
    } else if h < 1.0 / 2.0 {
        m2
    } else if h < 2.0 / 3.0 {
        m1 + (m2 - m1) * (2.0 / 3.0 - h) * 6.0
    } else {
        m1
    }
}

// ============================================================================
// Quantity_Color
// ============================================================================

/// A color stored as linear RGB triplet.
///
/// OCCT: `Quantity_Color`.
///
/// Internal storage is linear RGB in f32, matching OCCT's `NCollection_Vec3<float>`.
/// Constructors accept multiple color spaces; values are always converted to linear RGB.
fn color_epsilon() -> &'static std::sync::Mutex<f64> {
    static EPS: std::sync::LazyLock<std::sync::Mutex<f64>> =
        std::sync::LazyLock::new(|| std::sync::Mutex::new(1e-4));
    &EPS
}

#[derive(Debug, Clone)]
pub struct Color {
    rgb: [f32; 3],
}

impl Color {
    /// OCCT: default constructor → Quantity_NOC_YELLOW.
    pub fn new() -> Self {
        Self::from_name(NameOfColor::Yellow)
    }

    /// OCCT: `Quantity_Color(theName)` — from named color.
    pub fn from_name(name: NameOfColor) -> Self {
        Self {
            rgb: named_color_rgb(name),
        }
    }

    /// OCCT: `Quantity_Color(theC1, theC2, theC3, theType)`.
    pub fn from_values(c1: f64, c2: f64, c3: f64, color_type: TypeOfColor) -> Self {
        let rgb = match color_type {
            TypeOfColor::RGB => [c1 as f32, c2 as f32, c3 as f32],
            TypeOfColor::sRGB => srgb_to_linear_v3([c1 as f32, c2 as f32, c3 as f32]),
            TypeOfColor::HLS => {
                let srgb = hls_to_srgb([c1 as f32, c2 as f32, c3 as f32]);
                srgb_to_linear_v3(srgb)
            }
            TypeOfColor::CIELab => {
                // CIELab → XYZ → linear RGB
                lab_to_linear_rgb(c1, c2, c3)
            }
            TypeOfColor::CIELch => {
                let (c, h) = (c2, c3);
                let a = c * (h * std::f64::consts::PI / 180.0).cos();
                let b = c * (h * std::f64::consts::PI / 180.0).sin();
                lab_to_linear_rgb(c1, a, b)
            }
        };
        Self { rgb }
    }

    /// OCCT: `Quantity_Color(NCollection_Vec3<float>)` — from linear RGB vec directly.
    pub fn from_rgb_vec(rgb: [f32; 3]) -> Self {
        Self { rgb }
    }

    /// OCCT: `Red()`.
    pub fn red(&self) -> f64 {
        self.rgb[0] as f64
    }
    /// OCCT: `Green()`.
    pub fn green(&self) -> f64 {
        self.rgb[1] as f64
    }
    /// OCCT: `Blue()`.
    pub fn blue(&self) -> f64 {
        self.rgb[2] as f64
    }
    /// Access linear RGB as f32 slice.
    pub fn rgb(&self) -> [f32; 3] {
        self.rgb
    }

    /// OCCT: `Delta(theColor, DC, DI)` — percentage change of contrast and intensity.
    pub fn delta(&self, other: &Color) -> (f64, f64) {
        let hls_self = srgb_to_hls(linear_to_srgb_v3(self.rgb));
        let hls_other = srgb_to_hls(linear_to_srgb_v3(other.rgb));

        let dc = if hls_self[2].abs() > 1e-12 {
            ((hls_other[2] - hls_self[2]) / hls_self[2] * 100.0) as f64
        } else {
            0.0
        };
        let di = if hls_self[1].abs() > 1e-12 {
            ((hls_other[1] - hls_self[1]) / hls_self[1] * 100.0) as f64
        } else {
            0.0
        };
        (dc, di)
    }

    /// OCCT: `Distance(theColor)` — Euclidean distance in RGB space, range [0, sqrt(3)].
    pub fn distance(&self, other: &Color) -> f64 {
        let dr = self.rgb[0] as f64 - other.rgb[0] as f64;
        let dg = self.rgb[1] as f64 - other.rgb[1] as f64;
        let db = self.rgb[2] as f64 - other.rgb[2] as f64;
        (dr * dr + dg * dg + db * db).sqrt()
    }

    /// OCCT: `SquareDistance(theColor)`.
    pub fn square_distance(&self, other: &Color) -> f64 {
        let dr = self.rgb[0] as f64 - other.rgb[0] as f64;
        let dg = self.rgb[1] as f64 - other.rgb[1] as f64;
        let db = self.rgb[2] as f64 - other.rgb[2] as f64;
        dr * dr + dg * dg + db * db
    }

    /// OCCT: `IsEqual(theOther)`.
    pub fn is_equal(&self, other: &Color) -> bool {
        self.square_distance(other) <= Self::epsilon() * Self::epsilon()
    }

    /// OCCT: `IsDifferent(theOther)`.
    pub fn is_different(&self, other: &Color) -> bool {
        !self.is_equal(other)
    }

    /// OCCT: `Hue()` — returns hue in degrees [0, 360] or -1 for grayscale.
    pub fn hue(&self) -> f64 {
        srgb_to_hls(linear_to_srgb_v3(self.rgb))[0] as f64
    }

    /// OCCT: `Light()` — lightness component [0, 1].
    pub fn light(&self) -> f64 {
        srgb_to_hls(linear_to_srgb_v3(self.rgb))[1] as f64
    }

    /// OCCT: `Saturation()` — saturation component [0, 1].
    pub fn saturation(&self) -> f64 {
        srgb_to_hls(linear_to_srgb_v3(self.rgb))[2] as f64
    }

    /// OCCT: `ChangeIntensity(theDelta)` — delta as percentage of current lightness.
    pub fn change_intensity(&mut self, delta: f64) {
        let mut hls = srgb_to_hls(linear_to_srgb_v3(self.rgb));
        hls[1] = (hls[1] as f64 * (1.0 + delta / 100.0)).max(0.0).min(1.0) as f32;
        let srgb = hls_to_srgb(hls);
        self.rgb = srgb_to_linear_v3(srgb);
    }

    /// OCCT: `ChangeContrast(theDelta)` — delta as percentage of current saturation.
    pub fn change_contrast(&mut self, delta: f64) {
        let mut hls = srgb_to_hls(linear_to_srgb_v3(self.rgb));
        hls[2] = (hls[2] as f64 * (1.0 + delta / 100.0)).max(0.0).min(1.0) as f32;
        let srgb = hls_to_srgb(hls);
        self.rgb = srgb_to_linear_v3(srgb);
    }

    /// OCCT: `Name()` — find nearest named color.
    pub fn name(&self) -> NameOfColor {
        let mut best = NameOfColor::Black;
        let mut best_dist = f64::INFINITY;
        for nc in NAMED_COLORS.iter().copied() {
            let dist = self.square_distance(&Color {
                rgb: named_color_rgb(nc),
            });
            if dist < best_dist {
                best_dist = dist;
                best = nc;
            }
        }
        best
    }

    /// OCCT: `ColorFromName(theName, theColor)`.
    pub fn from_string_name(name: &str) -> Option<Self> {
        let upper = name.to_uppercase();
        for nc in NAMED_COLORS.iter().copied() {
            if nc.string_name() == upper {
                return Some(Self::from_name(nc));
            }
        }
        None
    }

    /// OCCT: `ColorFromHex(theHexColorString, theColor)`.
    /// Supports "#RGB", "#RRGGBB" format.
    pub fn from_hex(hex: &str) -> Option<Self> {
        let s = hex.trim_start_matches('#');
        let (r, g, b) = match s.len() {
            3 => {
                let r = u8::from_str_radix(&s[0..1], 16).ok()? * 17;
                let g = u8::from_str_radix(&s[1..2], 16).ok()? * 17;
                let b = u8::from_str_radix(&s[2..3], 16).ok()? * 17;
                (r, g, b)
            }
            6 => {
                let r = u8::from_str_radix(&s[0..2], 16).ok()?;
                let g = u8::from_str_radix(&s[2..4], 16).ok()?;
                let b = u8::from_str_radix(&s[4..6], 16).ok()?;
                (r, g, b)
            }
            _ => return None,
        };
        // Hex strings are sRGB → convert to linear
        let srgb = [r as f32 / 255.0, g as f32 / 255.0, b as f32 / 255.0];
        Some(Self {
            rgb: srgb_to_linear_v3(srgb),
        })
    }

    /// OCCT: `ColorToHex(theColor)`.
    pub fn to_hex(&self) -> String {
        let srgb = linear_to_srgb_v3(self.rgb);
        let r = (srgb[0] * 255.0 + 0.5) as u8;
        let g = (srgb[1] * 255.0 + 0.5) as u8;
        let b = (srgb[2] * 255.0 + 0.5) as u8;
        format!("#{:02X}{:02X}{:02X}", r, g, b)
    }

    /// OCCT: `Epsilon()`.
    pub fn epsilon() -> f64 {
        *color_epsilon().lock().unwrap()
    }

    /// OCCT: `SetEpsilon(theEpsilon)`.
    pub fn set_epsilon(eps: f64) {
        *color_epsilon().lock().unwrap() = eps;
    }

    /// OCCT: `StringName(theColor)` — name of a named color.
    pub fn string_name(name: NameOfColor) -> &'static str {
        name.string_name()
    }

    /// Create from 8-bit sRGB values.
    pub fn from_rgb8(r: u8, g: u8, b: u8) -> Self {
        let srgb = [r as f32 / 255.0, g as f32 / 255.0, b as f32 / 255.0];
        Self {
            rgb: srgb_to_linear_v3(srgb),
        }
    }

    /// Convert to ARGB int (OCCT Color2argb). Alpha = 0.
    pub fn to_argb(&self) -> i32 {
        let srgb = linear_to_srgb_v3(self.rgb);
        let r = (srgb[0] * 255.0 + 0.5) as i32;
        let g = (srgb[1] * 255.0 + 0.5) as i32;
        let b = (srgb[2] * 255.0 + 0.5) as i32;
        ((r & 0xff) << 16) | ((g & 0xff) << 8) | (b & 0xff)
    }

    /// Create from ARGB int (OCCT Argb2color). Alpha ignored.
    pub fn from_argb(argb: i32) -> Self {
        let r = ((argb >> 16) & 0xff) as u8;
        let g = ((argb >> 8) & 0xff) as u8;
        let b = (argb & 0xff) as u8;
        Self::from_rgb8(r, g, b)
    }
}

impl Default for Color {
    fn default() -> Self {
        Self::new()
    }
}

impl PartialEq for Color {
    fn eq(&self, other: &Self) -> bool {
        self.is_equal(other)
    }
}

impl Eq for Color {}

impl fmt::Display for Color {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Color(r={:.4}, g={:.4}, b={:.4})", self.red(), self.green(), self.blue())
    }
}

// ============================================================================
// ColorRGBA (OCCT Quantity_ColorRGBA)
// ============================================================================

/// Color with alpha component.
///
/// OCCT: `Quantity_ColorRGBA`.
#[derive(Debug, Clone, PartialEq)]
pub struct ColorRGBA {
    color: Color,
    alpha: f32,
}

impl ColorRGBA {
    /// Default: opaque yellow.
    pub fn new() -> Self {
        Self {
            color: Color::new(),
            alpha: 1.0,
        }
    }

    /// From Color + alpha.
    pub fn from_color(color: Color, alpha: f32) -> Self {
        Self { color, alpha }
    }

    /// From RGBA f32 values (linear RGB).
    pub fn from_rgba(r: f32, g: f32, b: f32, a: f32) -> Self {
        Self {
            color: Color::from_rgb_vec([r, g, b]),
            alpha: a,
        }
    }

    /// From sRGB 8-bit values + alpha.
    pub fn from_rgba8(r: u8, g: u8, b: u8, a: u8) -> Self {
        Self {
            color: Color::from_rgb8(r, g, b),
            alpha: a as f32 / 255.0,
        }
    }

    pub fn get_rgb(&self) -> &Color {
        &self.color
    }

    pub fn alpha(&self) -> f32 {
        self.alpha
    }

    pub fn set_alpha(&mut self, alpha: f32) {
        self.alpha = alpha;
    }

    pub fn is_different(&self, other: &ColorRGBA) -> bool {
        self.color.is_different(&other.color)
            || (self.alpha - other.alpha).abs() > Color::epsilon() as f32
    }
}

impl Default for ColorRGBA {
    fn default() -> Self {
        Self::new()
    }
}

// ============================================================================
// CIELab → linear RGB conversion (OCCT internals)
// ============================================================================

/// CIELab → XYZ → linear RGB (D65 white point).
fn lab_to_linear_rgb(l: f64, a: f64, b: f64) -> [f32; 3] {
    // Lab → XYZ (D65)
    let fy = (l + 16.0) / 116.0;
    let fx = a / 500.0 + fy;
    let fz = fy - b / 200.0;

    let x = lab_component(fx) * 0.95047; // D65 X
    let y = lab_component(fy);
    let z = lab_component(fz) * 1.08883; // D65 Z

    // XYZ → linear sRGB (D65)
    let r = x * 3.2406 + y * (-1.5372) + z * (-0.4986);
    let g = x * (-0.9689) + y * 1.8758 + z * 0.0415;
    let bv = x * 0.0557 + y * (-0.2040) + z * 1.0570;

    [r.clamp(0.0, 1.0) as f32, g.clamp(0.0, 1.0) as f32, bv.clamp(0.0, 1.0) as f32]
}

fn lab_component(t: f64) -> f64 {
    if t > 24.0 / 116.0 {
        t * t * t
    } else {
        108.0 / 841.0 * (t - 4.0 / 29.0)
    }
}

// All named colors in sorted order for Name() lookup.
const NAMED_COLORS: &[NameOfColor] = &[
    NameOfColor::AliceBlue,
    NameOfColor::AntiqueWhite,
    NameOfColor::Aquamarine,
    NameOfColor::Azure,
    NameOfColor::Beige,
    NameOfColor::Bisque,
    NameOfColor::Black,
    NameOfColor::BlanchedAlmond,
    NameOfColor::Blue,
    NameOfColor::BlueViolet,
    NameOfColor::Brown,
    NameOfColor::BurlyWood,
    NameOfColor::CadetBlue,
    NameOfColor::Chartreuse,
    NameOfColor::Chocolate,
    NameOfColor::Coral,
    NameOfColor::CornflowerBlue,
    NameOfColor::Cornsilk,
    NameOfColor::Cyan,
    NameOfColor::DarkGoldenrod,
    NameOfColor::DarkGreen,
    NameOfColor::DarkKhaki,
    NameOfColor::DarkOliveGreen,
    NameOfColor::DarkOrange,
    NameOfColor::DarkOrchid,
    NameOfColor::DarkRed,
    NameOfColor::DarkSalmon,
    NameOfColor::DarkSeaGreen,
    NameOfColor::DarkSlateBlue,
    NameOfColor::DarkSlateGray,
    NameOfColor::DarkTurquoise,
    NameOfColor::DarkViolet,
    NameOfColor::DeepPink,
    NameOfColor::DeepSkyBlue,
    NameOfColor::DimGray,
    NameOfColor::DodgerBlue,
    NameOfColor::Firebrick,
    NameOfColor::FloralWhite,
    NameOfColor::ForestGreen,
    NameOfColor::Fuchsia,
    NameOfColor::Gainsboro,
    NameOfColor::GhostWhite,
    NameOfColor::Gold,
    NameOfColor::Goldenrod,
    NameOfColor::Gray,
    NameOfColor::Green,
    NameOfColor::GreenYellow,
    NameOfColor::Honeydew,
    NameOfColor::HotPink,
    NameOfColor::IndianRed,
    NameOfColor::Indigo,
    NameOfColor::Ivory,
    NameOfColor::Khaki,
    NameOfColor::Lavender,
    NameOfColor::LavenderBlush,
    NameOfColor::LawnGreen,
    NameOfColor::LemonChiffon,
    NameOfColor::LightBlue,
    NameOfColor::LightCoral,
    NameOfColor::LightCyan,
    NameOfColor::LightGoldenrod,
    NameOfColor::LightGoldenrodYellow,
    NameOfColor::LightGray,
    NameOfColor::LightPink,
    NameOfColor::LightSalmon,
    NameOfColor::LightSeaGreen,
    NameOfColor::LightSkyBlue,
    NameOfColor::LightSlateBlue,
    NameOfColor::LightSlateGray,
    NameOfColor::LightSteelBlue,
    NameOfColor::LightYellow,
    NameOfColor::LimeGreen,
    NameOfColor::Linen,
    NameOfColor::Magenta,
    NameOfColor::Maroon,
    NameOfColor::MediumAquamarine,
    NameOfColor::MediumBlue,
    NameOfColor::MediumOrchid,
    NameOfColor::MediumPurple,
    NameOfColor::MediumSeaGreen,
    NameOfColor::MediumSlateBlue,
    NameOfColor::MediumSpringGreen,
    NameOfColor::MediumTurquoise,
    NameOfColor::MediumVioletRed,
    NameOfColor::MidnightBlue,
    NameOfColor::MintCream,
    NameOfColor::MistyRose,
    NameOfColor::Moccasin,
    NameOfColor::NavajoWhite,
    NameOfColor::NavyBlue,
    NameOfColor::OldLace,
    NameOfColor::Olive,
    NameOfColor::OliveDrab,
    NameOfColor::Orange,
    NameOfColor::OrangeRed,
    NameOfColor::Orchid,
    NameOfColor::PaleGoldenrod,
    NameOfColor::PaleGreen,
    NameOfColor::PaleTurquoise,
    NameOfColor::PaleVioletRed,
    NameOfColor::PapayaWhip,
    NameOfColor::PeachPuff,
    NameOfColor::Peru,
    NameOfColor::Pink,
    NameOfColor::Plum,
    NameOfColor::PowderBlue,
    NameOfColor::Purple,
    NameOfColor::Red,
    NameOfColor::RosyBrown,
    NameOfColor::RoyalBlue,
    NameOfColor::SaddleBrown,
    NameOfColor::Salmon,
    NameOfColor::SandyBrown,
    NameOfColor::SeaGreen,
    NameOfColor::Seashell,
    NameOfColor::Sienna,
    NameOfColor::Silver,
    NameOfColor::SkyBlue,
    NameOfColor::SlateBlue,
    NameOfColor::SlateGray,
    NameOfColor::Snow,
    NameOfColor::SpringGreen,
    NameOfColor::SteelBlue,
    NameOfColor::Tan,
    NameOfColor::Teal,
    NameOfColor::Thistle,
    NameOfColor::Tomato,
    NameOfColor::Turquoise,
    NameOfColor::Violet,
    NameOfColor::VioletRed,
    NameOfColor::Wheat,
    NameOfColor::White,
    NameOfColor::WhiteSmoke,
    NameOfColor::Yellow,
    NameOfColor::YellowGreen,
];

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_color() {
        let c = Color::new();
        assert!((c.red() - 1.0).abs() < 1e-5);
        assert!((c.green() - 1.0).abs() < 1e-5);
        assert_eq!(c.name(), NameOfColor::Yellow);
    }

    #[test]
    fn test_named_color() {
        let c = Color::from_name(NameOfColor::Red);
        assert!((c.red() - 1.0).abs() < 1e-5);
        assert!((c.green() - 0.0).abs() < 1e-5);
        assert!((c.blue() - 0.0).abs() < 1e-5);
    }

    #[test]
    fn test_rgb8() {
        let c = Color::from_rgb8(255, 0, 0);
        assert!((c.red() - 1.0).abs() < 1e-3);
    }

    #[test]
    fn test_hex_roundtrip() {
        let c = Color::from_hex("#FF0000").unwrap();
        let hex = c.to_hex();
        assert_eq!(hex, "#FF0000");
    }

    #[test]
    fn test_rgb_constructor() {
        let c = Color::from_values(0.5, 0.3, 0.1, TypeOfColor::RGB);
        assert!((c.red() - 0.5).abs() < 1e-5);
        assert!((c.green() - 0.3).abs() < 1e-5);
    }

    #[test]
    fn test_srgb_roundtrip() {
        let srgb_in = [0.8f32, 0.4, 0.2];
        let linear = srgb_to_linear_v3(srgb_in);
        let srgb_out = linear_to_srgb_v3(linear);
        for i in 0..3 {
            assert!((srgb_in[i] - srgb_out[i]).abs() < 1e-4);
        }
    }

    #[test]
    fn test_hls_roundtrip() {
        let c_orig = Color::from_values(0.7, 0.2, 0.5, TypeOfColor::RGB);
        let h = c_orig.hue();
        let l = c_orig.light();
        let s = c_orig.saturation();
        let c_restored = Color::from_values(h, l, s, TypeOfColor::HLS);
        let d = c_orig.distance(&c_restored);
        assert!(d < 1e-2, "HLS roundtrip distance too large: {}", d);
    }

    #[test]
    fn test_distance() {
        let black = Color::from_name(NameOfColor::Black);
        let white = Color::from_name(NameOfColor::White);
        assert!((black.distance(&white) - 1.732).abs() < 1e-3);

        let same = Color::from_name(NameOfColor::Red);
        assert!(black.is_different(&white));
        assert!(same.is_equal(&same));
    }

    #[test]
    fn test_hex_short() {
        let c = Color::from_hex("#F00").unwrap();
        assert!((c.red() - 1.0).abs() < 1e-3);
    }

    #[test]
    fn test_string_name() {
        let c = Color::from_string_name("RED").unwrap();
        assert!((c.red() - 1.0).abs() < 1e-3);
        assert!(Color::from_string_name("UNKNOWN").is_none());
    }

    #[test]
    fn test_intensity_contrast() {
        let mut c = Color::from_name(NameOfColor::Red);
        c.change_intensity(50.0);
        assert!(c.light() > 0.4); // should be lighter
    }

    #[test]
    fn test_argb() {
        let c = Color::from_name(NameOfColor::Red);
        let argb = c.to_argb();
        let c2 = Color::from_argb(argb);
        assert!(c.is_equal(&c2));
    }

    #[test]
    fn test_rgba() {
        let rgba = ColorRGBA::from_rgba8(255, 0, 0, 128);
        assert!((rgba.alpha() - 0.5).abs() < 0.01);
        assert!((rgba.get_rgb().red() - 1.0).abs() < 1e-3);
    }

    #[test]
    fn test_file_epsilon() {
        let old = Color::epsilon();
        Color::set_epsilon(1e-5);
        assert!((Color::epsilon() - 1e-5).abs() < 1e-12);
        Color::set_epsilon(old);
    }
}
