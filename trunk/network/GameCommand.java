package network;

import java.io.*;

import config.*;

/**
 * A <code>GameCommand</code> is a<code>Command</code> that is send by the
 * server to all the connected clients. This command contains the commands
 * that the clients have previousely send to the server.
 * This command is send when the server has recieved a <code>Command</code>
 * from all the players.
 */
public class GameCommand extends Command
{
    // Something
    public static final long serialVersionUID = 1;

    // The commands from all the clients
    private Command[] commands;

    /**
     * Creates a new <code>GameCommand</code> with an array of
     * the <code>Command</code>s from all the connected clients.
     *
     * @param commands An array containing the commands
     *                 from all the connected clients.
     */
    public GameCommand(Command[] commands)
    {
        // Sender is always server, so sender_id is 0
        super(Config.CMD_GAME, Config.SERVER_ID);
        this.commands = commands;
    }

    /**
     * Returns the <code>Commands</code> in this <code>GameCommand</code>.
     *
     * @param return The <code>Commands</code> in this command.
     */
    public Command[] getCommands()
    {
        return commands;
    }
}