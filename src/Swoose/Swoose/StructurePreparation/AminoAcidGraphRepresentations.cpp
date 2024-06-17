/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "AminoAcidGraphRepresentations.h"
#include <Molassembler/Graph.h>
#include <Molassembler/Molecule.h>
#include <Molassembler/Serialization.h>

namespace Scine {
namespace StructurePreparation {
namespace AminoAcids {

static const std::map<std::string, std::string> mapOfGraphStrings{
    {"ALA", "o2FhgqRhYQBhYwJhcqNhbIKBAYEDYmxygoEAgQFhc4KBAYEDYXMBpGFhAWFjAWFyo2Fsg4EAgQKBBGJscoOBAYECgQBhc4OBAoEEgQBhcw"
            "NhZ6JhRYSDAAEAgwECAIMBBACDAgMAYVqFBwYGCAZhdoMBAAA="},
    {"ARG",
     "o2Fhh6RhYQBhYwhhcqNhbIOBB4EJgQpibHKCgQCCAQJhc4KBB4IJCmFzAqRhYQBhYwJhcqNhbIKBAYEDYmxygoEAgQFhc4KBAYEDYXMBpGFhAGFjA"
     "WFyo2Fsg4EAgQKBBGJscoOBAoEBgQBhc4OBBIECgQBhcwOkYWEAYWMEYXKjYWyCgQGBBWJscoKBAYEAYXOCgQWBAWFzAaRhYQBhYwVhcqNhbIKBBI"
     "EGYmxygoEAgQFhc4KBBIEGYXMBpGFhAGFjBmFyo2FsgoEFgQdibHKCgQCBAWFzgoEFgQdhcwGkYWEAYWMHYXKjYWyCgQaBCGJscoKBAIEBYXOCgQa"
     "BCGFzAWFnomFFioMAAQCDAQIAgwEEAIMCAwCDBAUAgwUGAIMGBwCDBwgAgwgJAIMICgBhWosHBgYIBgYGBwYHB2F2gwEAAA=="},
    {"ASN", "o2FhhKRhYQBhYwVhcqNhbIOBBIEGgQdibHKDgQCBAoEBYXODgQSBB4EGYXMCpGFhAGFjBGFyo2FsgoEBgQVibHKCgQCBAWFzgoEBgQVhcw"
            "GkYWEAYWMBYXKjYWyDgQCBAoEEYmxyg4ECgQGBAGFzg4EEgQKBAGFzA6RhYQBhYwJhcqNhbIKBAYEDYmxygoEAgQFhc4KBAYEDYXMBYWei"
            "YUWHgwABAIMBAgCDAQQAgwIDAIMEBQCDBQYAgwUHAGFaiAcGBggGBggHYXaDAQAA"},
    {"ASP", "o2FhhKRhYQBhYwRhcqNhbIOBA4EFgQZibHKCgQCCAQJhc4KBA4IFBmFzAqRhYQBhYwNhcqNhbIKBAYEEYmxygoEAgQFhc4KBAYEEYXMBpG"
            "FhAGFjAWFyo2Fsg4EAgQKBA2JscoOBAoEBgQBhc4OBA4ECgQBhcwOkYWEAYWMCYXKjYWyCgQGBB2JscoKBAIEBYXOCgQGBB2FzAWFnomFF"
            "h4MAAQCDAQIAgwEDAIMCBwCDAwQAgwQFAIMEBgBhWogHBgYGBggICGF2gwEAAA=="},
    {"CYS", "o2Fhg6RhYQBhYwRhcqNhbIKBAYEFYmxygoEAgQFhc4KBAYEFYXMBpGFhAWFjAWFyo2Fsg4EAgQKBBGJscoOBAYECgQBhc4OBAoEEgQBhcw"
            "OkYWEAYWMCYXKjYWyCgQGBA2JscoKBAIEBYXOCgQGBA2FzAWFnomFFhYMAAQCDAQIAgwEEAIMCAwCDBAUAYVqGBwYGCAYQYXaDAQAA"},
    {"GLN",
     "o2FhhaRhYQBhYwVhcqNhbIOBBIEGgQdibHKDgQCBAoEBYXODgQSBB4EGYXMCpGFhAGFjBGFyo2FsgoEDgQVibHKCgQCBAWFzgoEDgQVhcwGkYWEAY"
     "WMDYXKjYWyCgQGBBGJscoKBAYEAYXOCgQSBAWFzAaRhYQBhYwFhcqNhbIOBAIECgQNibHKDgQKBAYEAYXODgQOBAoEAYXMDpGFhAGFjAmFyo2Fsgo"
     "EBgQhibHKCgQCBAWFzgoEBgQhhcwFhZ6JhRYiDAAEAgwECAIMBAwCDAggAgwMEAIMEBQCDBQYAgwUHAGFaiQcGBgYGBggHCGF2gwEAAA=="},
    {"GLU",
     "o2FhhaRhYQBhYwVhcqNhbIOBBIEGgQdibHKCgQCCAQJhc4KBBIIGB2FzAqRhYQBhYwRhcqNhbIKBA4EFYmxygoEAgQFhc4KBA4EFYXMBpGFhAGFjA"
     "2Fyo2FsgoEBgQRibHKCgQGBAGFzgoEEgQFhcwGkYWEAYWMBYXKjYWyDgQCBAoEDYmxyg4ECgQGBAGFzg4EDgQKBAGFzA6RhYQBhYwJhcqNhbIKBAY"
     "EIYmxygoEAgQFhc4KBAYEIYXMBYWeiYUWIgwABAIMBAgCDAQMAgwIIAIMDBACDBAUAgwUGAIMFBwBhWokHBgYGBgYICAhhdoMBAAA="},
    {"GLY", "o2FhgqRhYQBhYwJhcqNhbIKBAYEDYmxygoEAgQFhc4KBAYEDYXMBpGFhAGFjAWFyo2FsgoEAgQJibHKCgQGBAGFzgoECgQBhcwFhZ6JhRY"
            "ODAAEAgwECAIMCAwBhWoQHBgYIYXaDAQAA"},
    {"HIS", "pGFhiKRhYQBhYwhhcqRhbIKBBoEHY2xua4GiYXCCAAFjc2VxhQgGBAUHYmxygoEAgQFhc4KBBoEHYXMBpGFhAGFjB2FypGFsgoEFgQhjbG"
            "5rgaJhcIIAAWNzZXGFBwUEBghibHKCgQCBAWFzgoEFgQhhcwGkYWEAYWMCYXKjYWyCgQGBCWJscoKBAIEBYXOCgQGBCWFzAaRhYQBhYwFh"
            "cqNhbIOBAIECgQNibHKDgQKBAYEAYXODgQOBAoEAYXMDpGFhAGFjA2Fyo2FsgoEBgQRibHKCgQGBAGFzgoEEgQFhcwGkYWEAYWMEYXKkYW"
            "yDgQOBBYEGY2xua4GiYXCCAQJjc2VxhQQFBwgGYmxyg4EAgQKBAWFzg4EDgQaBBWFzAqRhYQBhYwVhcqRhbIKBBIEHY2xua4GiYXCCAAFj"
            "c2VxhQUEBggHYmxygoEAgQFhc4KBBIEHYXMBpGFhAGFjBmFypGFsgoEEgQhjbG5rgaJhcIIAAWNzZXGFBgQFBwhibHKCgQCBAWFzgoEEgQ"
            "hhcwFhYoWiYWEAYWWCBwiiYWEAYWWCBgiiYWEAYWWCBAWiYWEAYWWCBQeiYWEAYWWCBAZhZ6JhRYqDAAEAgwECAIMBAwCDAgkAgwMEAIME"
            "BQCDBAYAgwUHAIMGCACDBwgAYVqKBwYGBgYHBgYHCGF2gwEAAA=="},
    {"ILE", "o2FhhKRhYQBhYwRhcqNhbIKBA4EGYmxygoEAgQFhc4KBA4EGYXMBpGFhAGFjA2Fyo2Fsg4EBgQSBBWJscoOBAYEAgQJhc4OBBIEBgQVhcw"
            "OkYWEAYWMBYXKjYWyDgQCBAoEDYmxyg4ECgQGBAGFzg4EDgQKBAGFzA6RhYQBhYwJhcqNhbIKBAYEHYmxygoEAgQFhc4KBAYEHYXMBYWei"
            "YUWHgwABAIMBAgCDAQMAgwIHAIMDBACDAwUAgwQGAGFaiAcGBgYGBgYIYXaDAQAA"},
    {"LEU", "o2FhhKRhYQBhYwRhcqNhbIOBA4EFgQZibHKCgQCCAQJhc4KBA4IFBmFzA6RhYQBhYwNhcqNhbIKBAYEEYmxygoEBgQBhc4KBBIEBYXMBpG"
            "FhAGFjAWFyo2Fsg4EAgQKBA2JscoOBAoEBgQBhc4OBA4ECgQBhcwOkYWEAYWMCYXKjYWyCgQGBB2JscoKBAIEBYXOCgQGBB2FzAWFnomFF"
            "h4MAAQCDAQIAgwEDAIMCBwCDAwQAgwQFAIMEBgBhWogHBgYGBgYGCGF2gwEAAA=="},
    {"LYS", "o2FhhqRhYQBhYwZhcqNhbIKBBYEHYmxygoEAgQFhc4KBBYEHYXMBpGFhAGFjBWFyo2FsgoEEgQZibHKCgQCBAWFzgoEEgQZhcwGkYWEAYW"
            "MEYXKjYWyCgQOBBWJscoKBAIEBYXOCgQOBBWFzAaRhYQBhYwNhcqNhbIKBAYEEYmxygoEBgQBhc4KBBIEBYXMBpGFhAGFjAWFyo2Fsg4EA"
            "gQKBA2JscoOBAoEBgQBhc4OBA4ECgQBhcwOkYWEAYWMCYXKjYWyCgQGBCGJscoKBAIEBYXOCgQGBCGFzAWFnomFFiIMAAQCDAQIAgwEDAI"
            "MCCACDAwQAgwQFAIMFBgCDBgcAYVqJBwYGBgYGBgcIYXaDAQAA"},
    {"MET", "o2FhhaRhYQBhYwVhcqNhbIKBBIEGYmxygoEAgQFhc4KBBIEGYXMBpGFhAGFjBGFyo2FsgoEDgQVibHKCgQCBAWFzgoEDgQVhcwGkYWEAYW"
            "MDYXKjYWyCgQGBBGJscoKBAIEBYXOCgQGBBGFzAaRhYQBhYwFhcqNhbIOBAIECgQNibHKDgQKBAYEAYXODgQOBAoEAYXMDpGFhAGFjAmFy"
            "o2FsgoEBgQdibHKCgQCBAWFzgoEBgQdhcwFhZ6JhRYeDAAEAgwECAIMBAwCDAgcAgwMEAIMEBQCDBQYAYVqIBwYGBgYQBghhdoMBAAA="},
    {"PHE",
     "pGFhiaRhYQBhYwlhcqRhbIKBB4EIY2xua4GiYXCCAAFjc2VxhgkHBQQGCGJscoGCAAFhc4GCBwhhcwGkYWEAYWMIYXKkYWyCgQaBCWNsbmuBomFwg"
     "gABY3NlcYYIBgQFBwlibHKCgQCBAWFzgoEGgQlhcwGkYWEAYWMHYXKkYWyCgQWBCWNsbmuBomFwggABY3NlcYYHBQQGCAlibHKCgQCBAWFzgoEFgQ"
     "lhcwGkYWEAYWMCYXKjYWyCgQGBCmJscoKBAIEBYXOCgQGBCmFzAaRhYQBhYwFhcqNhbIOBAIECgQNibHKDgQKBAYEAYXODgQOBAoEAYXMDpGFhAGF"
     "jA2Fyo2FsgoEBgQRibHKCgQGBAGFzgoEEgQFhcwGkYWEAYWMEYXKkYWyDgQOBBYEGY2xua4GiYXCCAQJjc2VxhgQFBwkIBmJscoKCAQKBAGFzgoIF"
     "BoEDYXMCpGFhAGFjBWFypGFsgoEEgQdjbG5rgaJhcIIAAWNzZXGGBQQGCAkHYmxygoEAgQFhc4KBBIEHYXMBpGFhAGFjBmFypGFsgoEEgQhjbG5rg"
     "aJhcIIAAWNzZXGGBgQFBwkIYmxygoEAgQFhc4KBBIEIYXMBYWKGomFhAGFlgggJomFhAGFlggcJomFhAGFlggYIomFhAGFlggQFomFhAGFlggUHom"
     "FhAGFlggQGYWeiYUWLgwABAIMBAgCDAQMAgwIKAIMDBACDBAUAgwQGAIMFBwCDBggAgwcJAIMICQBhWosHBgYGBgYGBgYGCGF2gwEAAA=="},
    {"PRO", "o2FhhqRhYQBhYwVhcqRhbIKBAIEEY2xua4GiYXCCAAFjc2VxhQUAAQMEYmxygoEBgQBhc4KBBIEAYXMBpGFhAGFjBGFypGFsgoEDgQVjbG"
            "5rgaJhcIIAAWNzZXGFBAMBAAVibHKCgQCBAWFzgoEDgQVhcwGkYWEAYWMDYXKkYWyCgQGBBGNsbmuBomFwggABY3NlcYUDAQAFBGJscoKB"
            "AYEAYXOCgQSBAWFzAaRhYQBhYwJhcqNhbIKBAYEGYmxygoEAgQFhc4KBAYEGYXMBpGFhAGFjAGFypGFsgoEBgQVjbG5rgaJhcIIAAWNzZX"
            "GFAAEDBAVibHKCgQCBAWFzgoEBgQVhcwGkYWEAYWMBYXKkYWyDgQCBAoEDY2xua4GiYXCCAAJjc2VxhQEABQQDYmxyg4ECgQGBAGFzg4ED"
            "gQKBAGFzA2FnomFFh4MAAQCDAAUAgwECAIMBAwCDAgYAgwMEAIMEBQBhWocHBgYGBgYIYXaDAQAA"},
    {"SER", "o2Fhg6RhYQBhYwNhcqNhbIKBAYEEYmxygoEAgQFhc4KBAYEEYXMBpGFhAGFjAWFyo2Fsg4EAgQKBA2JscoKCAQKBAGFzgoICA4EAYXMDpG"
            "FhAGFjAmFyo2FsgoEBgQVibHKCgQCBAWFzgoEBgQVhcwFhZ6JhRYWDAAEAgwECAIMBAwCDAgUAgwMEAGFahgcGBgYICGF2gwEAAA=="},
    {"THR",
     "o2Fhg6RhYQBhYwNhcqNhbIOBAYEEgQVibHKDgQCBAoEBYXODgQGBBYEEYXMDpGFhAGFjAWFyo2Fsg4EAgQKBA2JscoOBAoEBgQBhc4OBA4ECgQBhc"
     "wOkYWEAYWMCYXKjYWyCgQGBBmJscoKBAIEBYXOCgQGBBmFzAWFnomFFhoMAAQCDAQIAgwEDAIMCBgCDAwQAgwMFAGFahwcGBgYIBghhdoMBAAA="},
    {"TRP", "pGFhjKRhYQBhYwxhcqRhbIKBCoELY2xua4GiYXCCAAFjc2VxhgwKCAYJC2JscoKBAYEAYXOCgQuBCmFzAaRhYQBhYwthcqRhbIKBCYEMY2"
            "xua4GiYXCCAAFjc2VxhgsJBggKDGJscoKBAIEBYXOCgQmBDGFzAaRhYQBhYwphcqRhbIKBCIEMY2xua4GiYXCCAAFjc2VxhgoIBgkLDGJs"
            "coKBAYEAYXOCgQyBCGFzAaRhYQBhYwlhcqRhbIKBBoELY2xua4GiYXCCAAFjc2VxhgkGCAoMC2JscoKBAIEBYXOCgQaBC2FzAaRhYQBhYw"
            "hhcqRhbIOBBoEHgQpjbG5rgqJhcIIAAWNzZXGFCAYEBQeiYXCCAAJjc2VxhggGCQsMCmJscoOBAIECgQFhc4OBBoEKgQdhcwKkYWEAYWMH"
            "YXKkYWyCgQWBCGNsbmuBomFwggABY3NlcYUHBQQGCGJscoKBAYEAYXOCgQiBBWFzAaRhYQBhYwJhcqNhbIKBAYENYmxygoEAgQFhc4KBAY"
            "ENYXMBpGFhAGFjAWFyo2Fsg4EAgQKBA2JscoOBAoEBgQBhc4OBA4ECgQBhcwOkYWEAYWMDYXKjYWyCgQGBBGJscoKBAYEAYXOCgQSBAWFz"
            "AaRhYQBhYwRhcqRhbIOBA4EFgQZjbG5rgaJhcIIBAmNzZXGFBAUHCAZibHKDgQKBAIEBYXODgQaBA4EFYXMCpGFhAGFjBWFypGFsgoEEgQ"
            "djbG5rgaJhcIIAAWNzZXGFBQQGCAdibHKCgQCBAWFzgoEEgQdhcwGkYWEAYWMGYXKkYWyDgQSBCIEJY2xua4KiYXCCAAFjc2VxhQYEBQcI"
            "omFwggECY3NlcYYGCAoMCwlibHKDgQCBAoEBYXODgQSBCYEIYXMCYWKKomFhAGFlggsMomFhAGFlggkLomFhAGFlggoMomFhAGFlgggKom"
            "FhAGFlggQGomFhAGFlggUHomFhAGFlggQFomFhAGFlggYIomFhAGFlggcIomFhAGFlggYJYWeiYUWPgwABAIMBAgCDAQMAgwINAIMDBACD"
            "BAUAgwQGAIMFBwCDBggAgwYJAIMHCACDCAoAgwkLAIMKDACDCwwAYVqOBwYGBgYGBgcGBgYGBghhdoMBAAA="},
    {"VAL",
     "o2Fhg6RhYQBhYwNhcqNhbIOBAYEEgQVibHKCgQCCAQJhc4KBAYIEBWFzA6RhYQBhYwFhcqNhbIOBAIECgQNibHKDgQKBAYEAYXODgQOBAoEAYXMDp"
     "GFhAGFjAmFyo2FsgoEBgQZibHKCgQCBAWFzgoEBgQZhcwFhZ6JhRYaDAAEAgwECAIMBAwCDAgYAgwMEAIMDBQBhWocHBgYGBgYIYXaDAQAA"},
    {"TYR", "pGFhiaRhYQBhYwlhcqRhbIOBB4EIgQpjbG5rgaJhcIIAAWNzZXGGCQcFBAYIYmxygoIAAYECYXOCggcIgQphcwKkYWEAYWMIYXKkYWyCgQ"
            "aBCWNsbmuBomFwggABY3NlcYYIBgQFBwlibHKCgQCBAWFzgoEGgQlhcwGkYWEAYWMHYXKkYWyCgQWBCWNsbmuBomFwggABY3NlcYYHBQQG"
            "CAlibHKCgQCBAWFzgoEFgQlhcwGkYWEAYWMCYXKjYWyCgQGBC2JscoKBAIEBYXOCgQGBC2FzAaRhYQBhYwFhcqNhbIOBAIECgQNibHKDgQ"
            "KBAYEAYXODgQOBAoEAYXMDpGFhAGFjA2Fyo2FsgoEBgQRibHKCgQGBAGFzgoEEgQFhcwGkYWEAYWMEYXKkYWyDgQOBBYEGY2xua4GiYXCC"
            "AQJjc2VxhgQFBwkIBmJscoKCAQKBAGFzgoIFBoEDYXMCpGFhAGFjBWFypGFsgoEEgQdjbG5rgaJhcIIAAWNzZXGGBQQGCAkHYmxygoEAgQ"
            "Fhc4KBBIEHYXMBpGFhAGFjBmFypGFsgoEEgQhjbG5rgaJhcIIAAWNzZXGGBgQFBwkIYmxygoEAgQFhc4KBBIEIYXMBYWKGomFhAGFlgggJ"
            "omFhAGFlggcJomFhAGFlggYIomFhAGFlggQFomFhAGFlggUHomFhAGFlggQGYWeiYUWMgwABAIMBAgCDAQMAgwILAIMDBACDBAUAgwQGAI"
            "MFBwCDBggAgwcJAIMICQCDCQoAYVqMBwYGBgYGBgYGBggIYXaDAQAA"},
    // Dicysteine
    {"DCY", "o2FhiKRhYQBhYwthcqNhbIKBBYEKYmxygoEBgQBhc4KBCoEFYXMBpGFhAGFjCmFyo2FsgoEHgQtibHKCgQCBAWFzgoEHgQthcwGkYWEAYW"
            "MCYXKjYWyCgQGBA2JscoKBAIEBYXOCgQGBA2FzAaRhYQFhYwFhcqNhbIOBAIECgQRibHKDgQGBAoEAYXODgQKBBIEAYXMDpGFhAGFjCGFy"
            "o2FsgoEHgQlibHKCgQCBAWFzgoEHgQlhcwGkYWEAYWMEYXKjYWyCgQGBBWJscoKBAIEBYXOCgQGBBWFzAaRhYQBhYwVhcqNhbIKBBIELYm"
            "xygoEAgQFhc4KBBIELYXMBpGFhAWFjB2Fyo2Fsg4EGgQiBCmJscoOBAYECgQBhc4OBCIEKgQZhcwNhZ6JhRYuDAAEAgwECAIMBBACDAgMA"
            "gwQFAIMFCwCDBgcAgwcIAIMHCgCDCAkAgwoLAGFajAcGBggGEAcGBggGEGF2gwEAAA=="},
    // Selenomethionine
    {"MSE", "o2FhhaRhYQBhYwdhcqNhbIKBBYEGYmxygoEAgQFhc4KBBYEGYXMBpGFhAGFjBWFyo2FsgoEEgQdibHKCgQCBAWFzgoEEgQdhcwGkYWEAYW"
            "MEYXKjYWyCgQGBBWJscoKBAIEBYXOCgQGBBWFzAaRhYQBhYwFhcqNhbIOBAIECgQRibHKDgQKBAYEAYXODgQSBAoEAYXMDpGFhAGFjAmFy"
            "o2FsgoEBgQNibHKCgQCBAWFzgoEBgQNhcwFhZ6JhRYeDAAEAgwECAIMBBACDAgMAgwQFAIMFBwCDBgcAYVqIBwYGCAYGBhgiYXaDAQAA"},
    // Selenocysteine
    {"SEC",
     "o2Fhg6RhYQBhYwRhcqNhbIKBAYEFYmxygoEAgQFhc4KBAYEFYXMBpGFhAWFjAWFyo2Fsg4EAgQKBBGJscoOBAYECgQBhc4OBAoEEgQBhcwOkYWEAY"
     "WMCYXKjYWyCgQGBA2JscoKBAIEBYXOCgQGBA2FzAWFnomFFhYMAAQCDAQIAgwEEAIMCAwCDBAUAYVqGBwYGCAYYImF2gwEAAA=="},
    // Pyrrolysine
    {"PYL", "o2FhjaRhYQBhYw9hcqRhbIKBAoEMY2xua4GiYXCCAAFjc2VxhQ8CCAoMYmxygoEAgQFhc4KBAoEMYXMBpGFhAGFjDmFyo2Fsg4ECgQ2BEG"
            "JscoOBAIEBgQJhc4OBAoENgRBhcwKkYWEAYWMNYXKjYWyCgQuBDmJscoKBAIEBYXOCgQuBDmFzAaRhYQBhYwxhcqRhbIKBCoEPY2xua4Gi"
            "YXCCAAFjc2VxhQwKCAIPYmxygoEAgQFhc4KBCoEPYXMBpGFhAGFjC2Fyo2FsgoEJgQ1ibHKCgQCBAWFzgoEJgQ1hcwGkYWEAYWMKYXKkYW"
            "yCgQiBDGNsbmuBomFwggABY3NlcYUKCAIPDGJscoKBAIEBYXOCgQiBDGFzAaRhYQBhYwlhcqNhbIKBB4ELYmxygoEAgQFhc4KBB4ELYXMB"
            "pGFhAWFjAmFypGFsg4EIgQ6BD2NsbmuBomFwggACY3NlcYUCCAoMD2JscoOBAIEBgQJhc4OBCIEOgQ9hcwOkYWEAYWMBYXKjYWyDgQCBA4"
            "EFYmxyg4ECgQGBAGFzg4EFgQOBAGFzA6RhYQFhYwhhcqRhbIOBAoEGgQpjbG5rgaJhcIIAAmNzZXGFCAIPDApibHKDgQKBAIEBYXODgQqB"
            "AoEGYXMDpGFhAGFjA2Fyo2FsgoEBgQRibHKCgQCBAWFzgoEBgQRhcwGkYWEAYWMFYXKjYWyCgQGBB2JscoKBAYEAYXOCgQeBAWFzAaRhYQ"
            "BhYwdhcqNhbIKBBYEJYmxygoEAgQFhc4KBBYEJYXMBYWeiYUWRgwABAIMBAwCDAQUAgwIIAIMCDgCDAg8AgwMEAIMFBwCDBggAgwcJAIMI"
            "CgCDCQsAgwoMAIMLDQCDDA8Agw0OAIMOEABhWpEHBgYGCAYGBgYGBgYGBwYHCGF2gwEAAA=="},
};

// Transform the graph's string representation to a Molassembler::Graph object
std::map<std::string, Molassembler::Graph> getAminoAcidGraphs() {
  std::map<std::string, Molassembler::Graph> mapOfAminoAcidGraphs;
  for (auto const& pair : mapOfGraphStrings) {
    const std::string serializationString = pair.second;
    auto binary = Molassembler::JsonSerialization::base64Decode(serializationString);
    Molassembler::JsonSerialization serialization(binary, Molassembler::JsonSerialization::BinaryFormat::CBOR);
    Molassembler::Molecule aminoAcidMolecule = serialization;
    auto aminoAcidGraph = aminoAcidMolecule.graph();
    mapOfAminoAcidGraphs.insert({pair.first, aminoAcidGraph});
  }
  return mapOfAminoAcidGraphs;
}

// Alanine
const std::vector<std::string> ALATypes = {"N", "CA", "C", "O", "CB"};

// Arginine
const std::vector<std::string> ARGTypes = {"N", "CA", "C", "O", "CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"};

// Asparagine
const std::vector<std::string> ASNTypes = {"N", "CA", "C", "O", "CB", "CG", "OD1", "ND2"};

// Aspartic acid
const std::vector<std::string> ASPTypes = {"N", "CA", "C", "CB", "CG", "OD1", "OD2", "O"};

// Cysteine
const std::vector<std::string> CYSTypes = {"N", "CA", "C", "O", "CB", "SG"};

// Glutamine
const std::vector<std::string> GLNTypes = {"N", "CA", "C", "CB", "CG", "CD", "OE1", "NE2", "O"};

// Glutamic acid
const std::vector<std::string> GLUTypes = {"N", "CA", "C", "CB", "CG", "CD", "OE1", "OE2", "O"};

// Glycine
const std::vector<std::string> GLYTypes = {"N", "CA", "C", "O"};

// Histidine
const std::vector<std::string> HISTypes = {"N", "CA", "C", "CB", "CG", "ND1", "CD2", "CE1", "NE2", "O"};

// Isoleucine
const std::vector<std::string> ILETypes = {"N", "CA", "C", "CB", "CG1", "CG2", "CD1", "O"};

// Leucine
const std::vector<std::string> LEUTypes = {"N", "CA", "C", "CB", "CG", "CD1", "CD2", "O"};

// Lysine
const std::vector<std::string> LYSTypes = {"N", "CA", "C", "CB", "CG", "CD", "CE", "NZ", "O"};

// Methionine
const std::vector<std::string> METTypes = {"N", "CA", "C", "CB", "CG", "SD", "CE", "O"};

// Phenylalanine
const std::vector<std::string> PHETypes = {"N", "CA", "C", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "O"};

// Proline
const std::vector<std::string> PROTypes = {"N", "CA", "C", "CB", "CG", "CD", "O"};

// Serine
const std::vector<std::string> SERTypes = {"N", "CA", "C", "CB", "OG", "O"};

// Threonin
const std::vector<std::string> THRTypes = {"N", "CA", "C", "CB", "OG1", "CG2", "O"};

// Tryptophan
const std::vector<std::string> TRPTypes = {"N",   "CA",  "C",   "CB",  "CG",  "CD1", "CD2",
                                           "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2", "O"};

// Tyrosine
const std::vector<std::string> TYRTypes = {"N", "CA", "C", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH", "O"};

// Valine
const std::vector<std::string> VALTypes = {"N", "CA", "C", "CB", "CG1", "CG2", "O"};

// Dicysteine
const std::vector<std::string> DCYTypes = {"N", "CA", "C", "O", "CB", "SG", "N", "CA", "C", "O", "CB", "SG"};

// Selenomethionine
const std::vector<std::string> MSETypes = {"N", "CA", "C", "O", "CB", "CG", "SE", "CE"};

// Selenocysteine
const std::vector<std::string> SECTypes = {"N", "CA", "C", "O", "CB", "SE"};

// Pyrrolysine
const std::vector<std::string> PYLTypes = {"N",  "CA",  "CA2", "C",   "O",  "CB", "CB2", "CG", "CG2",
                                           "CD", "CD2", "CE",  "CE2", "NZ", "C2", "N2",  "O2"};

std::vector<std::string> getResidueTypes(const std::string& residueName) {
  if (residueName == "ALA")
    return ALATypes;
  else if (residueName == "ARG")
    return ARGTypes;
  else if (residueName == "ASN")
    return ASNTypes;
  else if (residueName == "ASP")
    return ASPTypes;
  else if (residueName == "CYS")
    return CYSTypes;
  else if (residueName == "GLN")
    return GLNTypes;
  else if (residueName == "GLU")
    return GLUTypes;
  else if (residueName == "GLY")
    return GLYTypes;
  else if (residueName == "HIS")
    return HISTypes;
  else if (residueName == "ILE")
    return ILETypes;
  else if (residueName == "LEU")
    return LEUTypes;
  else if (residueName == "LYS")
    return LYSTypes;
  else if (residueName == "MET")
    return METTypes;
  else if (residueName == "PHE")
    return PHETypes;
  else if (residueName == "PRO")
    return PROTypes;
  else if (residueName == "SER")
    return SERTypes;
  else if (residueName == "THR")
    return THRTypes;
  else if (residueName == "TRP")
    return TRPTypes;
  else if (residueName == "TYR")
    return TYRTypes;
  else if (residueName == "VAL")
    return VALTypes;
  else if (residueName == "DCY")
    return DCYTypes;
  else if (residueName == "MSE")
    return MSETypes;
  else if (residueName == "SEC")
    return SECTypes;
  else if (residueName == "PYL")
    return PYLTypes;
  else
    throw std::runtime_error("Unsupported Amino Acid");
}

} // namespace AminoAcids
} // namespace StructurePreparation
} // namespace Scine
