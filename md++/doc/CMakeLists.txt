

if(DOXYGEN_FOUND)

    # Check if Python3 is installed
    find_package(Python3 REQUIRED COMPONENTS Interpreter)

    # Configure Doxygen
    set(DOXYGEN_INPUT_FILE ${CMAKE_CURRENT_SOURCE_DIR}/doxygen.conf.in)
    set(DOXYGEN_OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/doxygen.conf)

    configure_file(${DOXYGEN_INPUT_FILE} ${DOXYGEN_OUTPUT_FILE} @ONLY)

    # Generate snippets
    set(SNIPPETS_DIR ${CMAKE_CURRENT_BINARY_DIR}/snippets)
    file(MAKE_DIRECTORY ${SNIPPETS_DIR})

    add_custom_command(
        OUTPUT ${SNIPPETS_DIR}/snippets.cc
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/generate_snippets.py
                ${SNIPPETS_DIR}/snippets.cc
                ${CMAKE_CURRENT_SOURCE_DIR}/../src/io/parameter/in_parameter.cc
                ${CMAKE_CURRENT_SOURCE_DIR}/../src/io/topology/in_*.cc
        DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/generate_snippets.py
        COMMENT "Generating documentation snippets"
    )

    # Generate documentation
    add_custom_target(doc
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUTPUT_FILE}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM
    )

    # Ensure the `doc` target depends on the snippets being generated
    add_custom_target(generate_snippets DEPENDS ${SNIPPETS_DIR}/snippets.cc)
    add_dependencies(doc generate_snippets)

    # Copy figures to the build directory
    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/figures
         DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/doc)

    # Install the generated documentation
    install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc/html
            DESTINATION ${CMAKE_INSTALL_PREFIX}/doc)

    # Install additional documentation files
    install(FILES
        doxygen.conf.in
        install.doxy
        install.html
        main.doxy
        DESTINATION ${CMAKE_INSTALL_PREFIX}/doc
    )

    # Optional: Generate PDF documentation
    add_custom_target(doc_pdf
        COMMAND pdflatex -output-directory ${CMAKE_CURRENT_BINARY_DIR}/doc ${CMAKE_CURRENT_SOURCE_DIR}/gromosXX.tex
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc
        COMMENT "Generating PDF documentation"
        VERBATIM
    )
else()
    message(WARNING "Doxygen not found. Documentation will not be generated.")
endif()