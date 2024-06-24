function mesh = generateRectangleMesh(options)
    arguments 
        options.width = 100
        options.height = 20
        options.nElementsX = 20
        options.nElementsZ = 4
        options.filename = 'tmp.mat'
        options.pythonExecutable = '/usr/local/bin/python3'
    end

    pyenv('Version',options.pythonExecutable);
    pyrunfile(sprintf('rectangle.py --filename %s --width %d --height %d --nElementsZ %i --nElementsX %i', options.filename, options.width, options.height, options.nElementsZ, options.nElementsX));
    
    mesh = load(options.filename);

    mesh = SFMesh(mesh.nodes, mesh.elements, mesh.boundaries);

end