####### Install rules
include(GNUInstallDirs)

# Install executable
install(TARGETS MRCDFT
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} # bin
    COMPONENT EXE
)

# Install documentation
install(FILES
    README.md
    DESTINATION ${CMAKE_INSTALL_DOCDIR} # sahre/doc/MRCDFT
    COMPONENT DOCS
)