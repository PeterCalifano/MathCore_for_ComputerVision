# MathCore is header-only: install this module's public headers without adding
# compiled sources to the parent library inventory.
file(GLOB installable_headers "*.h")

file(RELATIVE_PATH relative_path_from_root
     "${PROJECT_SOURCE_DIR}"
     "${CMAKE_CURRENT_SOURCE_DIR}")
string(REPLACE "src/" "" relative_path_without_src "${relative_path_from_root}")

install(FILES ${installable_headers}
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${project_name}/${relative_path_without_src}")
