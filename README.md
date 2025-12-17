# Maze Generator
This is a simple maze generator that allows you to generate any maze you'd like using c++.

Create a maze by running the executable and passing the name of the desired mask as the first argument. For example: *g++ a.exe bird*

A few templates are provided in the repository, and a user can create their own templates that use the following rules:
- The template is plain text in a .txt file
- The template must "draw" a shape in a monospaced font using the '1' character.
- Each row of characters must be seperated with an explicit new line (not by automatic text wrapping).
- The shape must be contiguous. This means that the shape must be in "one piece", and that you can travel between any two characters by travelling across a series of adjacent characters.

As of now, new masks must be created by hand. A possible future project could be to create masks from images directly, or to accept mask files with less strict guidelines.
