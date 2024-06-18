import fs from 'fs';
import path from 'path';
    
const base_path = '..';

function scanDirectory(directory, result = []) {
    const files = fs.readdirSync(directory);
    console.log(`Scanning ${directory}`);

    files.forEach((file) => {
        const filePath = path.join(directory, file);
        const stats = fs.statSync(filePath);

        if (stats.isDirectory()) {
            scanDirectory(filePath, result);
        } else {
            result.push({filename: file, path: directory, source: filePath.replace(publicFolder, base_path)});
        }
    });

    return result;
}

const publicFolder = 'public';
const fileTree = scanDirectory(publicFolder);
const jsonOutput = JSON.stringify(fileTree, null, 2);

fs.writeFileSync(path.join('public', 'preloaded_files.json'), jsonOutput);