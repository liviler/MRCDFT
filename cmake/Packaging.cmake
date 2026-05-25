####### CPack configuration
# Package basic info
set(CPACK_PACKAGE_NAME "MRCDFT")
set(CPACK_PACKAGE_VERSION "1.0.0")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Multi-reference Covariant Density Functional Theory code")

# Output directory
set(CPACK_OUTPUT_FILE_PREFIX "${CMAKE_BINARY_DIR}/package")

# Install prefix inside package
# set(CPACK_PACKAGING_INSTALL_PREFIX "/opt/mrcdft")

# Files to include in package
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/LICENSE.txt")
set(CPACK_RESOURCE_FILE_README "${CMAKE_SOURCE_DIR}/README.md")



# Generator selection (Windows/Linux auto handled)
if(WIN32)
    set(CPACK_GENERATOR "ZIP;NSIS;WIX")

    set(CPACK_PACKAGE_INSTALL_DIRECTORY "MRCDFT")

    ## Only for NSIS
    # Enable NSIS installer option to modify PATH
    set(CPACK_NSIS_MODIFY_PATH ON)

    ## Only for WIX
    set(CPACK_WIX_UPGRADE_GUID
        "8A3F4E3A-7B1A-4D2B-9C91-123456789ABC"
    )
    set(WIX_PATCH_DIR "${CMAKE_CURRENT_SOURCE_DIR}/cmake/wix")
    set(CPACK_WIX_MRCDFT_ID "CM_CP_EXE.bin.MRCDFT.exe")
    configure_file(
        "${WIX_PATCH_DIR}/wix_patch.xml.in"
        "${WIX_PATCH_DIR}/wix_patch.xml"
        @ONLY
      )
    set(CPACK_WIX_PATCH_FILE "${WIX_PATCH_DIR}/wix_patch.xml")

else()
    set(CPACK_GENERATOR "TGZ")
endif()

# Components
set(CPACK_COMPONENTS_ALL EXE DOCS)

# Tell CPack what belongs to which component
set(CPACK_COMPONENT_EXE_DISPLAY_NAME "MRCDFT Core")
set(CPACK_COMPONENT_DOCS_DISPLAY_NAME "Documentation")
include(CPack)