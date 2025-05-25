Coding Conventions
==================

Indentation and Formatting
--------------------------

- Use 4 spaces for indentation
- No tabs
- Maximum line length is 80 characters
- Never produce trailing spaces on a line (unless we have a very special
  case like a test where the command produces trailing spaces)
- Opening braces on the same line as control statements
- Except in function definitions, there the opening brace is on a new line
- Closing braces on their own line
- Space after keywords like ``if``, ``for``, ``while``
- Space around operators (``=``, ``+``, ``-``, etc.)
- Space between function name and opening parenthesis
- Space after commas in parameter lists
- We generally have a space before opening parentheses or brackets, except
  for a *second* bracket (in multi-dimensional dereferences) e.g.::

    my_function ();
    my_array [0][5] = 5;

Naming Conventions
------------------

- Function names use PascalCase with PGA prefix for public API (e.g.,
  ``PGACreate``)
- Internal functions may use ``camelCase`` or ``snake_case`` but new function
  should use ``snake_case``
- Constants use ``ALL_CAPS`` with ``PGA_`` prefix (e.g., ``PGA_TRUE``)
- Variables use ``camelCase`` or ``snake_case``, newer code should use
  ``snake_case``
- Struct members use ``snake_case``

Function Structure
------------------

- Functions have a brief description in comments that can be extracted
  with doxygen
- Parameters and return values are documented
- Debug statements at entry and exit of functions using ``PGADebugEntered``
  and ``PGADebugExited`` for existing function, new functions will not use this
- Error checking with ``PGAError`` or ``PGAErrorPrintf``

Macros and Constants
--------------------

- Macros are ``ALL_CAPS``
- Constants are defined with ``#define``
- Enums are used for related constants

Comments
--------

- Function headers use a specific format with brief, description,
  parameters, return
- Internal comments use ``/* */`` style
- Temporary debug code may be commented with ``/* */``

Conditional Compilation
-----------------------

- Use ``#if !defined(DOXYGEN_SHOULD_SKIP_THIS)`` for sections to be skipped
  by documentation generator
- Use ``#ifdef``, ``#ifndef`` for platform-specific code

Error Handling
--------------

- ``PGAError`` and ``PGAErrorPrintf`` are equivalent, for new code
  ``PGAErrorPrintf`` is preferred
- Use ``PGAErrorPrintf`` with the PGA_FATAL flag (second parameter) for
  fatal errors
- Check parameters at function entry for functions in public API

Memory Management
-----------------

- Dynamic arrays on the stack with dynamic parameters are declared with
  ``DECLARE_DYNARRAY`` macro
- Pointer arithmetic used for array access

Version Control
---------------

- Commit messages start with a capitalized (50 chars or less) summary
  line, avoid ending the summary line with a period
- Then a longer description may follow after a blank line
- No line should be longer than 72 characters
- Commit messages do not use a type like "feat" or "docs" as mandated in
  other git conventions
- Use the present tense ("Add feature" not "Added feature")
- Use the imperative mood ("Move cursor to..." not "Moves cursor to...")
